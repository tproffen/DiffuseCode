MODULE structur

USE errlist_mod 
!
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
USE discus_config_mod 
USE discus_allocate_appl_mod
USE crystal_mod 
USE chem_mod 
USE diffuse_mod
USE molecule_mod 
USE prop_para_mod 
USE read_internal_mod
USE discus_save_mod 
USE spcgr_apply
USE spcgr_mod 
USE stack_rese_mod
USE update_cr_dim_mod
use wyckoff_mod
!      USE interface_def
!
USE doact_mod 
USE do_wait_mod
USE build_name_mod
USE learn_mod 
USE lib_echo
USE lib_errlist_func
USE lib_help
USE lib_length
USE lib_macro_func
USE get_params_mod
USE class_macro_internal
USE precision_mod
USE prompt_mod 
USE take_param_mod
USE str_comp_mod
USE sup_mod
USE support_mod
!
IMPLICIT none 
!                                                                       
INTEGER , PARAMETER :: MAXW = 11
!                                                                       
CHARACTER(LEN=PREC_STRING)        :: line, zeile, cpara (maxw) 
CHARACTER(LEN=PREC_STRING)        :: strucfile 
CHARACTER(LEN=LEN(prompt)) :: orig_prompt
CHARACTER(LEN=5)           :: befehl 
INTEGER, DIMENSION(MAXW)   :: lpara
INTEGER          :: lp, length 
INTEGER          :: lstr, i, j, k
INTEGER          :: ianz, lbef
LOGICAL          :: lout 
REAL(KIND=PREC_DP)   , DIMENSION(maxw) ::  werte!, wwerte
INTEGER          ::   occupancy= 0          ! Apply occupancy upon read cell   ?
LOGICAL          :: l_identical= .FALSE.    ! Are atoms allowed to be identical?
LOGICAL          :: l_site     = .FALSE.    ! Treat atoms on different sites as different types?
REAL             :: r_identical = 1.0E-5
INTEGER, PARAMETER :: NOPTIONAL = 5
INTEGER, PARAMETER :: O_SETTING = 4
INTEGER, PARAMETER :: O_SITE    = 5
CHARACTER(LEN=   9)       , DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=PREC_STRING), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: NCALC = 1 ! Number of values to calculate 
!
!
!                                                                       
!
DATA oname  / 'radius' , 'occupancy', 'identical', 'setting', 'site'   /
DATA loname /  9       ,  9         ,  6         ,  7       ,  4       /
!
opara  =  (/ '1.0E-5'  , 'clear '   , 'none  '   , 'abc   ' , 'equal ' /)   ! Always provide fresh default values
lopara =  (/  6        ,  6         ,  6         ,  6       ,  6       /)
owerte =  (/  1.0E-5   ,  0.0       ,  0.0       ,  0.0     ,  0.0     /)
!
!
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
orig_prompt = prompt
prompt = prompt (1:len_str (prompt) ) //'/read' 
!
9000  CONTINUE
!
   CALL get_cmd (line, length, befehl, lbef, zeile, lp, prompt) 
!                                                                       
   IF(ier_num.ne.0) RETURN 
   IF(line (1:1)  == ' '.or.line (1:1)  == '#' .or.   & 
      line == char(13) .or. line(1:1) == '!'  ) GOTO 9000
!                                                                       
!------ execute a macro file                                            
!                                                                       
   IF (line (1:1) .eq.'@') THEN 
      IF (length.ge.2) THEN 
         line(1:length-1) = line(2:length)
         line(length:length) = ' '
         length = length - 1
         CALL file_kdo(line, length)
      ELSE 
         ier_num = - 13 
         ier_typ = ER_MAC 
      ENDIF 
!                                                                       
!------ Echo a string, just for interactive check in a macro 'echo'     
!                                                                       
      ELSEIF (str_comp (befehl, 'echo', 2, lbef, 4) ) THEN 
         CALL echo (zeile, lp) 
!                                                                       
!     execute command                                                   
!     help                                                              
!                                                                       
      ELSEIF (str_comp (befehl, 'help', 2, lbef, 4) .or.str_comp (befehl&
     &, '?   ', 1, lbef, 4) ) THEN                                      
         IF (zeile.eq.' '.or.zeile.eq.char (13) ) THEN 
            zeile = 'commands' 
            lp = lp + 8 
         ENDIF 
         IF (str_comp (zeile, 'errors', 2, lp, 6) ) THEN 
            lp = lp + 7 
            CALL do_hel ('discus '//zeile, lp) 
         ELSE 
            lp = lp + 12 
            CALL do_hel ('discus read '//zeile, lp) 
         ENDIF 
!                                                                       
!------  -----waiting for user input                                    
!                                                                       
      ELSEIF (str_comp (befehl, 'wait', 3, lbef, 4) ) THEN 
         CALL do_input (zeile, lp) 
!                                                                       
      ELSE 
!                                                                       
!     --all other commands                                              
!                                                                       
!     --get command parameters                                          
!                                                                       
         CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
         IF (ier_num.ne.0) THEN 
            GOTO 8888              ! Jump to handle error messages, amd macro conditions
         ENDIF 
         IF (ianz.ge.1) THEN 
!                                                                       
!     --Build file name                                                 
!                                                                       
            CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
            IF (ier_num.ne.0) THEN 
               GOTO 8888              ! Jump to handle error messages, amd macro conditions
            ENDIF 
         ENDIF 
!
!      --Get optional parameters
!
         opara  =  (/ '1.0E-5'  , 'clear '   , 'none  '   , 'abc   ' , 'equal ' /)   ! Always provide fresh default values
         lopara =  (/  6        ,  6         ,  6         ,  6       ,  6       /)
         owerte =  (/  1.0E-5   ,  0.0       ,  0.0       ,  0.0     ,  0.0     /)
         CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                           oname, loname, opara, lopara, lpresent, owerte)
         IF(ier_num/=0) GOTO 8888              ! Jump to handle error messages, amd macro conditions
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
         IF (str_comp (befehl, 'cell',  1, lbef, 4) .or. &
             str_comp (befehl, 'lcell', 1, lbef, 5)      ) THEN                                    
            l_identical = str_comp (opara(3), 'tolerate', 5, lopara(3), 8)
            r_identical = owerte(1)
            IF(opara(2)=='clear')  THEN
               occupancy = 0
            ELSEIF(opara(2)=='apply')  THEN
               occupancy = 1
            ELSEIF(opara(2)=='keep')  THEN
               occupancy = -1
            ELSE
               ier_num = -13
               ier_typ = ER_COMM
               RETURN
            ENDIF
            l_site = opara(O_SITE) == 'differ'
            CALL do_readcell(befehl,lbef,ianz, maxw, cpara, lpara, &
                             l_identical, r_identical, occupancy,  &
                             l_site)
!                                                                       
!     Free style editing of a structure 'free'                          
!                                                                       
         ELSEIF (str_comp (befehl, 'free', 1, lbef, 4) ) THEN 
!           CALL do_readfree(befehl,lbef,ianz, maxw, cpara, lpara)
            CALL do_readfree(            ianz, maxw, cpara, lpara, &
                             opara(O_SETTING), lopara(O_SETTING))
!                                                                       
!     read an old structure 'stru'                                      
!                       Ä                                                
         ELSEIF (str_comp (befehl, 'stru', 1, lbef, 4) ) THEN 
            IF (ianz.eq.1) THEN 
               CALL rese_cr 
               sav_r_ncell = .false. 
               strucfile = cpara (1)(1:lpara(1))
               l_site = opara(O_SITE) == 'differ'
               CALL do_readstru(strucfile, l_site)
               IF(ier_num /= 0) THEN
                  IF(ier_msg(3) == ' ') THEN
                     ier_msg(3) = strucfile(MAX(1,LEN_TRIM(strucfile)-LEN(ier_msg)):LEN_TRIM(strucfile))
                  ENDIF 
               ENDIF 
!                                                                       
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
            IF (ier_num.eq.0) THEN 
!                                                                       
!     ------reset microdomain status                                    
!                                                                       
               CALL do_stack_rese 
!              Flag that no Fourier has been calculated yet
               four_last = FOUR_NN
            ENDIF 
         ELSEIF (str_comp (befehl, 'exit', 1, lbef, 4) ) THEN 
            GOTO 9999 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
            GOTO 9999 
         ENDIF 
8888     CONTINUE    ! Target for errors, in order to handle these properly
         IF (ier_num.eq.0) THEN 
            IF(cr_syst==4) THEN
               WRITE (output_io, 1000) cr_spcgr, cr_spcgrno , cr_set, cr_spcgr_set
            ELSE
               WRITE (output_io, 2000) cr_spcgr, cr_spcgrno
            ENDIF
!.......calculate metric and reciprocal metric tensor,reciprocal lattice
!       constants and permutation tensors                               
            CALL setup_lattice (cr_a0, cr_ar, cr_eps, cr_gten, cr_reps, &
            cr_rten, cr_win, cr_wrez, cr_v, cr_vr, lout, cr_gmat,       &
            cr_fmat, cr_cartesian,                                      &
              cr_tran_g, cr_tran_gi, cr_tran_f, cr_tran_fi)
            IF (.not. (str_comp (befehl, 'cell',  1, lbef, 4) .or.      &
                       str_comp (befehl, 'lcell', 1, lbef, 5)     ) ) THEN
               CALL get_symmetry_matrices 
            ENDIF
         ELSE 
            CALL errlist 
            IF (ier_sta.ne.ER_S_LIVE) THEN 
               IF (lmakro .OR. lmakro_error) THEN  ! Error within macro or termination errror
                  IF(sprompt /= prompt ) THEN
                     ier_num = -10
                     ier_typ = ER_COMM
                     ier_msg(1) = ' Error occured in read menu'
                     prompt_status = PROMPT_ON 
                     prompt = orig_prompt
                     RETURN
                  ELSE
                     CALL macro_close 
                     prompt_status = PROMPT_ON 
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
            CALL no_error 
         ENDIF 
      ENDIF 
!                                                                       
 9999 CONTINUE 
!
DO i=1,3
  IF(cr_icc(i)==1) chem_period(i) = .FALSE.
ENDDO
!
prompt = orig_prompt
!                                                                       
 1000 FORMAT    (1x,a16,i5, ' Setting:',a3,2x, a16) 
 2000 FORMAT    (1x,a16,i5)
!
END SUBROUTINE read_struc                     
!
!*******************************************************************************
!
SUBROUTINE do_readcell(befehl,lbef,ianz, maxw, cpara, lpara, l_identical, &
                       r_identical, occupancy, l_site)
!
USE discus_allocate_appl_mod
USE chem_mod 
USE crystal_mod
USE diffuse_mod
USE molecule_mod
USE prop_para_mod
USE read_internal_mod
USE stack_rese_mod
USE update_cr_dim_mod
!
USE ber_params_mod
USE errlist_mod
USE precision_mod
USE str_comp_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*),                  INTENT(IN) :: befehl
INTEGER         ,                  INTENT(IN) :: lbef
INTEGER         ,                  INTENT(IN) :: ianz
INTEGER         ,                  INTENT(IN) :: MAXW
CHARACTER(LEN=*), DIMENSION(MAXW), INTENT(INOUT) :: cpara
INTEGER         , DIMENSION(MAXW), INTENT(INOUT) :: lpara
LOGICAL         ,                  INTENT(IN) :: l_identical
REAL            ,                  INTENT(IN) :: r_identical
INTEGER         ,                  INTENT(IN) :: occupancy
LOGICAL         ,                  INTENT(IN) :: l_site
!
CHARACTER(LEN=MAX(PREC_STRING,LEN(cpara))) :: strucfile
CHARACTER(LEN=MAX(PREC_STRING,LEN(cpara))) :: outfile
LOGICAL             :: need_alloc
INTEGER             :: i, j, k, l, n
INTEGER             :: iatom
INTEGER             :: ce_natoms
INTEGER             :: n_atom
INTEGER             :: n_gene
INTEGER             :: n_symm
INTEGER             :: n_mole
INTEGER             :: n_type
INTEGER             :: ncells
INTEGER, DIMENSION(3) :: local_icc
INTEGER, DIMENSION(3), PARAMETER :: one = (/ 1, 1, 1/)
REAL(KIND=PREC_DP)   , DIMENSION(MAXW) :: werte
REAL                :: r
!
!
IF (ianz.ge.1) THEN 
   cr_newtype = str_comp (befehl, 'cell', 1, lbef, 4) 
   CALL rese_cr 
   strucfile = cpara(1) (1:lpara(1))
   IF (ier_num.eq.0) THEN 
!                                                                       
!     --------if necessary get crystal size                             
!                                                                       
      IF (ianz.gt.1) THEN 
         cpara (1) = '0.0' 
         lpara (1) = 3 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.eq.0) THEN 
            DO i = 1, ianz - 1 
               cr_icc (i) = nint (werte (i + 1) ) 
            ENDDO 
         ENDIF 
      ENDIF 
      local_icc(:) = cr_icc(:)
      IF (ier_num.eq.0) THEN 
         internalcell:        IF ( str_comp(strucfile(1:8),'internal',8,8,8)) THEN
            CALL readcell_internal(strucfile)
         ELSE internalcell
!                        CALL test_file ( strucfile, natoms, nscats, n_mole, n_type, &
!                                         n_atom, -1 , .false.)
!                        IF (ier_num /= 0) THEN
!                           RETURN
!                        ENDIF
            CALL import_test(0, strucfile, outfile)
            IF(ier_num == 0) THEN
               strucfile = outfile
               CALL readcell (strucfile, l_identical, r_identical) 
               cr_icc(:) = local_icc(:)   ! Restore cr_icc in case molecules were read
            ENDIF
         ENDIF internalcell
!
         IF (ier_num.eq.0) THEN 
!                                                                       
!     ----------check whether total number of atoms fits into available 
!     ----------space                                                   
!                                                                       
            iatom = cr_icc (1) * cr_icc (2) * cr_icc (3) * cr_natoms
            IF (iatom.gt.nmax) THEN 
               CALL alloc_crystal ( MAXSCAT, INT(iatom * 1.1))
               IF (ier_num < 0 ) THEN
                  ier_msg(3) = strucfile(MAX(1,LEN_TRIM(strucfile)-LEN(ier_msg)):LEN_TRIM(strucfile))
                  RETURN                 ! Jump to handle error messages, amd macro conditions
!                 GOTO 8888              ! Jump to handle error messages, amd macro conditions
               ENDIF
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
            cr_ncreal  = 0   ! Non void atoms in unit cell
            DO n=1,cr_natoms
               IF(cr_at_lis(cr_iscat(n))/='VOID') cr_ncreal = cr_ncreal + 1
            ENDDO
            IF(l_site) CALL differ_site(cr_nscat, cr_ncatoms, ONE)
            cr_natoms = 0 
            DO k = 1, cr_icc (3) 
               DO j = 1, cr_icc (2) 
                  DO i = 1, cr_icc (1) 
                     DO n = 1, ce_natoms 
                        cr_natoms = cr_natoms + 1 
                        cr_iscat (cr_natoms) = cr_iscat (n) 
                        cr_pos (1, cr_natoms) = cr_pos (1, n) + REAL( i - 1)
                        cr_pos (2, cr_natoms) = cr_pos (2, n) + REAL( j - 1)
                        cr_pos (3, cr_natoms) = cr_pos (3, n) + REAL( k - 1)
                        cr_mole (cr_natoms) = cr_mole (n)
                        cr_prop (cr_natoms) = cr_prop (n)
                        cr_surf(:,cr_natoms) = cr_surf(:,n)
                        cr_magn(:,cr_natoms) = cr_magn(:,n)
                     ENDDO 
                  ENDDO 
               ENDDO 
            ENDDO 
!
!     ---------- Apply occupancy
!
            IF(occupancy == 0) THEN     ! Clear occupancies on read cell
               cr_occ(:) = 1.0
            ELSEIF(occupancy==1) THEN   ! Apply occupancies
               DO i=1, cr_natoms
                  IF(cr_occ(cr_iscat(i))<1.0) THEN
                     CALL RANDOM_NUMBER(r)
                     IF(r > cr_occ(cr_iscat(i))) THEN
                        cr_iscat(i) = 0
                        cr_prop (i)  = ibclr (cr_prop (i),  PROP_NORMAL)
                     ENDIF
                  ENDIF
               ENDDO
               cr_occ(:) = 1.0
            ELSEIF(occupancy==-1) THEN   ! Ignore occupancies but keep values
               CONTINUE
            ENDIF
!                                                                       
!     ----------Update crystal dimensions                               
!                                                                       
            CALL update_cr_dim 
            ncells = cr_icc (1) * cr_icc (2)* cr_icc (3)
!                                                                       
!     ----------If molecules were read                                  
!                                                                       
            IF (mole_num_mole.gt.0) THEN 
               need_alloc = .false.
               n_gene = MAX( 1, MOLE_MAX_GENE)
               n_symm = MAX( 1, MOLE_MAX_SYMM)
               n_mole =         MOLE_MAX_MOLE
               n_type =         MOLE_MAX_TYPE
               n_atom =         MOLE_MAX_ATOM
               IF (mole_num_mole* ncells                              >= MOLE_MAX_MOLE ) THEN
                  n_mole = mole_num_mole* ncells + 20
                  need_alloc = .true.
               ENDIF
               IF ((mole_off(mole_num_mole)+mole_len(mole_num_mole))*ncells >= MOLE_MAX_ATOM ) THEN
                  n_atom = (mole_off(mole_num_mole)+mole_len(mole_num_mole))*ncells + 200
                  need_alloc = .true.
               ENDIF
               IF ( need_alloc ) THEN
                  call alloc_molecule(n_gene, n_symm, n_mole, n_type, n_atom)
               ENDIF
               IF (mole_num_mole * cr_icc (1) * cr_icc (2)  &
                                 * cr_icc (3) .le.MOLE_MAX_MOLE) THEN         
                  mole_num_atom = mole_off (mole_num_mole)  &
                                + mole_len (mole_num_mole)                
                  l             = mole_num_mole 
                  mole_num_unit = mole_num_mole 
                  DO i = 2,cr_icc(1)*cr_icc(2)*cr_icc(3)                                
                     DO j = 1, mole_num_mole 
                        l = l + 1 
                        mole_len (l) = mole_len (j) 
                        mole_off (l) = mole_off (l -              &
                        mole_num_mole) + mole_num_atom            
                        mole_type (l) = mole_type (j) 
                        mole_char (l) = mole_char (j) 
                        mole_dens (l) = mole_dens (j) 
!                       mole_biso (l) = mole_biso (j) 
                        DO k = 1, mole_len (j) 
                           mole_cont(mole_off(l) + k) =&
                           mole_cont(mole_off(j) + k) + (i - 1) * ce_natoms   
                           iatom          = mole_cont (mole_off (l) + k)
                           cr_prop(iatom) = ibset(cr_prop(iatom),PROP_MOLECULE)
                           cr_mole(iatom) = l
                        ENDDO 
                     ENDDO 
                  ENDDO 
                  mole_num_mole = l 
                  mole_num_atom = mole_off (mole_num_mole) +&
                                  mole_len (mole_num_mole)                
               ELSE 
                  ier_num = - 65 
                  ier_typ = ER_APPL 
                  RETURN                 ! Jump to handle error messages, amd macro conditions
!                 GOTO 8888              ! Jump to handle error messages, amd macro conditions
               ENDIF 
            ENDIF 
!                                                                       
!     ----------Define initial crystal size in fractional coordinates   
!               cr_dim0(:,1) is often used as the coordinate of the lower left
!               unit cell. Distances are THEN calculated relative to this
!               unit cell to obtain relative unit cell numbers. If a large molecule
!               sticks out of the unit cell, although its center is within the 
!               unit cell, the offset was calculated wrong. cr_dim(:,2) is hardly used. 
!               To reflect the intention of cr_dim0(:,1) it is now calculated from cr_icc.
            DO l = 1, 3 
!              cr_dim0 (l, 1) = REAL(nint (cr_dim (l, 1) ) ) 
!              cr_dim0 (l, 2) = REAL(nint (cr_dim (l, 2) ) ) 
               IF(MOD(cr_icc(l),2)==0) THEN
                              
!                             cr_dim0 (l, 1) = FLOAT(-(cr_icc(l)-1)/2)
!                             cr_dim0 (l, 2) = FLOAT((cr_icc(l)+1)/2)+1
                  cr_dim0(l,1) = FLOAT(-cr_icc(l)/2)
                  cr_dim0(l,2) = cr_dim0(l,1) + cr_icc(l) - 1
               ELSE
!                 cr_dim0 (l, 1) = FLOAT(-(cr_icc(l)  )/2)
!                 cr_dim0 (l, 2) = FLOAT((cr_icc(l)  )/2)+1
                  cr_dim0(l,1) = FLOAT(-(cr_icc(l)+1)/2 + 1)
                  cr_dim0(l,2) = cr_dim0(l,1) + cr_icc(l) - 1
               ENDIF 
            ENDDO 
         ENDIF 
      ENDIF 
   ENDIF 
ENDIF 
IF (ier_num.eq.0) THEN 
!                                                                       
!     ------reset microdomain status                                    
!                                                                       
   CALL do_stack_rese 
!  Flag that no Fourier has been calculated yet
   four_last = FOUR_NN
ELSE
   IF(ier_msg(3) == ' ') THEN
      ier_msg(3) = strucfile(MAX(1,LEN_TRIM(strucfile)-LEN(ier_msg)):LEN_TRIM(strucfile))
   ENDIF 
ENDIF 
!
chem_purge = .FALSE.    ! No purge was done, period boundary is OK
chem_period(:) = .TRUE.
chem_quick     = .TRUE.
!
END SUBROUTINE do_readcell
!
!*******************************************************************************
!
SUBROUTINE do_readstru(strucfile, l_site)
!
! Do the full job for a 'read stru ' command
!
USE crystal_mod 
USE chem_mod 
USE molecule_mod 
USE read_internal_mod
!
USE errlist_mod
USE precision_mod
USE str_comp_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: strucfile
LOGICAL         , INTENT(IN)    :: l_site     ! Differ atoms on sites?
!
CHARACTER(LEN=MAX(PREC_STRING,LEN(strucfile))) :: outfile 
!
CALL rese_cr
internals:     IF ( str_comp(strucfile(1:8),'internal',8,8,8)) THEN
   CALL readstru_internal(strucfile) !, NMAX, MAXSCAT, MOLE_MAX_MOLE, &
!                       MOLE_MAX_TYPE, MOLE_MAX_ATOM )
   IF(ier_num/=0) RETURN
ELSE internals
   CALL import_test(1, strucfile, outfile)
   IF(ier_num == 0) THEN
      strucfile = outfile
   ELSE
      RETURN
   ENDIF
   CALL do_readstru_disk(strucfile, l_site)
!
   IF (ier_num /= 0) THEN
      RETURN                 ! Jump to handle error messages, amd macro conditions
   ENDIF
!
ENDIF internals
!
chem_purge = .FALSE.    ! No purge was done, period boundary is OK
CALL test_mole_gap
!
IF(chem_quick .AND. l_site) THEN
   CALL differ_site(cr_nscat, cr_ncatoms, cr_icc )
ENDIF
!
END SUBROUTINE do_readstru
!
!*******************************************************************************
!
SUBROUTINE do_readstru_disk(strucfile, l_site)
!
! Do the readstru part for a file on disk
!
USE crystal_mod 
USE chem_mod 
USE discus_allocate_appl_mod
USE discus_save_mod 
USE molecule_mod 
USE read_internal_mod
!
USE errlist_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: strucfile
LOGICAL         , INTENT(IN)    :: l_site     ! Differ atoms on sites?
!
INTEGER             :: natoms
INTEGER             :: nscats
INTEGER             :: n_mole
INTEGER             :: n_type
INTEGER             :: n_atom
INTEGER             :: i,l
integer, dimension(3) :: n_cells
!
!
CALL test_file(strucfile, natoms, nscats, n_mole, n_type, &
               n_atom, n_cells, -1 , .false.)
IF (ier_num /= 0) THEN
   RETURN                 ! Jump to handle error messages, amd macro conditions
ENDIF
!
IF(natoms > NMAX .or. nscats > MAXSCAT) THEN
   natoms = MAX(INT(natoms * 1.1), natoms + 10,NMAX)
   nscats = MAX(INT(nscats * 1.1), nscats + 2, MAXSCAT)
   CALL alloc_crystal (nscats, natoms)
   IF ( ier_num /= 0 ) THEN
      RETURN                 ! Jump to handle error messages, amd macro conditions
   ENDIF
ENDIF
IF(n_mole>MOLE_MAX_MOLE .or. n_type>MOLE_MAX_TYPE .or.   &
   n_atom>MOLE_MAX_ATOM                          ) THEN
   n_mole = MAX(n_mole +20 ,MOLE_MAX_MOLE)
   n_type = MAX(n_type +10 ,MOLE_MAX_TYPE)
   n_atom = MAX(n_atom +200,MOLE_MAX_ATOM)
   CALL alloc_molecule(1, 1,n_mole,n_type,n_atom)
   IF ( ier_num /= 0 )  THEN
      RETURN                 ! Jump to handle error messages, amd macro conditions
   ENDIF
ENDIF
!
CALL readstru(NMAX, MAXSCAT, strucfile, cr_name,        &
              cr_spcgr, cr_set, cr_a0, cr_win, cr_natoms, cr_nscat, cr_dw, cr_occ,    &
              cr_at_lis, cr_pos, cr_mole, cr_surf, cr_magn, cr_iscat, cr_prop, cr_dim, cr_magnetic, as_natoms, &
              as_at_lis, as_dw, as_pos, as_iscat, as_prop, sav_ncell,  &
              sav_r_ncell, sav_ncatoms, spcgr_ianz, spcgr_para)        
IF(ier_num.ne.0) THEN 
   RETURN                 ! Jump to handle error messages, amd macro conditions
ENDIF 
mole_num_atom = mole_off (mole_num_mole)  &  !Update number of atoms in molecules
                + mole_len (mole_num_mole)                
!                                                                       
!     ------Define initial crystal size in fractional coordinates       
!                                                                       
DO l = 1, 3 
   cr_dim0 (l, 1) = REAL(nint (cr_dim (l, 1) ) ) 
   cr_dim0 (l, 2) = REAL(nint (cr_dim (l, 2) ) ) 
ENDDO 
!                                                                       
!     ------The crystal size was read from the structure file           
!                                                                       
IF (sav_r_ncell) THEN 
   DO i = 1, 3 
      cr_icc (i) = sav_ncell (i) 
   ENDDO 
   cr_ncatoms = sav_ncatoms 
   cr_ncreal  = sav_ncatoms 
ELSE 
!                                                                       
!     ------Define initial crystal size in number of unit cells         
!                                                                       
   DO i = 1, 3 
      cr_icc(i) = MAX(1,INT(cr_dim(i,2) - cr_dim(i,1) + 1. ) )                                                
   ENDDO 
!                                                                       
!     ------Define (average) number of atoms per unit cell              
!                                                                       
   cr_ncatoms = MAX(1,cr_natoms / (cr_icc (1) * cr_icc (2)     &
                                 * cr_icc (3) ))
ENDIF 
!
IF(cr_natoms /= cr_icc(1)*cr_icc(2)*cr_icc(3)*cr_ncatoms) THEN
   chem_period(:) = .false.
   chem_quick     = .false.
ELSE
   chem_period(:) = .TRUE.
   chem_quick     = .TRUE.
ENDIF
!
END SUBROUTINE do_readstru_disk
!
!********************************************************************** 
!
SUBROUTINE differ_site(n_types, n_ncatoms, icc)
!
USE crystal_mod
USE discus_allocate_appl_mod
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: n_types      ! Number of atom types
INTEGER, INTENT(IN) :: n_ncatoms    ! Number of atomRs per unit cell
INTEGER, DIMENSION(3), INTENT(IN) :: icc
INTEGER :: i,isite, k, iat, iscat
INTEGER :: l
INTEGER :: new_scat
INTEGER :: n_scat, n_max
CHARACTER(LEN=4), DIMENSION(0:n_types, n_ncatoms) :: names
INTEGER         , DIMENSION(0:n_types, n_ncatoms) :: types
INTEGER         , DIMENSION(0:n_types, n_ncatoms) :: old
INTEGER         , DIMENSION(           n_ncatoms) :: n_on_site
REAL            , DIMENSION(0:n_types, n_ncatoms) :: dbw
!
types = -1
n_on_site = 0
new_scat = 0
!
! Determine number of atom types on each site
!
cells: DO i=1, icc(1)*icc(2)*icc(3)
   sites: DO isite=1, n_ncatoms
      iat = (i-1)*n_ncatoms + isite
      iscat = cr_iscat(iat)
      DO k=0, n_on_site(isite)
         IF(cr_at_lis(iscat)==names(k,isite) .AND. cr_dw(iscat)==dbw(k,isite)) THEN
            CYCLE sites
         ENDIF
      ENDDO
      n_on_site(isite) = n_on_site(isite) + 1
      names(n_on_site(isite), isite) = cr_at_lis(iscat)
      dbw  (n_on_site(isite), isite) = cr_dw    (iscat)
      old  (n_on_site(isite), isite) = iscat
   ENDDO sites
ENDDO cells
!
! assign new atom type numbers
!
iscat = 0
DO l=1,cr_nscat                              ! Loop over all old scattering types
maxsites: DO i=1,n_types
   sites_b: DO isite=1, cr_ncatoms
      IF(i>n_on_site(isite)) CYCLE sites_b
      IF(l==old(i, isite)) THEN
      iscat = iscat + 1
      types(i, isite) = iscat
      ENDIF
   ENDDO sites_b
ENDDO maxsites
ENDDO
new_scat = iscat    ! New maximum atom types
!
! Assign new atom type to each atom
!
cells_c: DO i=1, icc(1)*icc(2)*icc(3)
   sites_c: DO isite=1, n_ncatoms
      iat = (i-1)*n_ncatoms + isite
      iscat = cr_iscat(iat)
      DO k=0, n_on_site(isite)
         IF(iscat==old(k,isite)) THEN
            cr_iscat(iat) = types(k,isite)
            CYCLE sites_c
         ENDIF
      ENDDO
   ENDDO sites_c
ENDDO cells_c
!
! Make new list of atom names and DBW's
!
IF(new_scat>MAXSCAT) THEN
   n_scat = new_scat
   n_max  = NMAX
   CALL alloc_crystal(n_scat, n_max)
ENDIF
!
DO isite=1,n_ncatoms
   DO i=1,n_on_site(isite)
      iscat = types(i, isite)
      cr_at_lis(iscat) = names(i, isite)
      cr_dw    (iscat) = dbw  (i, isite)
   ENDDO
ENDDO
cr_nscat= new_scat
!
END SUBROUTINE differ_site
!********************************************************************** 
!
SUBROUTINE do_readfree(ianz, maxw, cpara, lpara, c_set, l_set)
!
USE chem_mod
USE crystal_mod
USE diffuse_mod
USE spcgr_mod
USE spcgr_apply
USE stack_rese_mod
USE ber_params_mod
USE get_params_mod
USE lib_errlist_func
USE precision_mod
!
!CHARACTER(LEN=*),                  INTENT(IN) :: befehl
!INTEGER         ,                  INTENT(IN) :: lbef
INTEGER         ,                  INTENT(INOUT) :: ianz
INTEGER         ,                  INTENT(IN   ) :: MAXW
CHARACTER(LEN=*), DIMENSION(MAXW), INTENT(INOUT) :: cpara
INTEGER         , DIMENSION(MAXW), INTENT(INOUT) :: lpara
CHARACTER(LEN=*),                  INTENT(IN   ) :: c_set
INTEGER         ,                  INTENT(IN   ) :: l_set
!
INTEGER  :: i
INTEGER  :: iianz
!
REAL(KIND=PREC_DP), DIMENSION(MAXW) :: werte
REAL(KIND=PREC_DP), DIMENSION(MAXW) :: wwerte
!
CALL rese_cr 
cr_name = 'freely created structure' 
cr_spcgr (1:1)  = 'P' 
cr_spcgr (2:2)  = '1' 
cr_spcgr (3:16) = '              ' 
cr_set          = c_set(1:MIN(3,l_set))
cr_spcgrno = 1 
cr_syst = 1 
spcgr_para = 1
CALL get_symmetry_matrices    ! Define matrices for default P1

CALL no_error                 ! Ignore errors for the moment
IF (ianz==0) THEN 
   cr_a0 (:) = 1.0 
   cr_win(:) = 90.0 
ELSEIF (ianz==6 .OR. ianz==7 .OR. ianz==8) THEN 
   iianz = 6
   CALL ber_params (iianz, cpara, lpara, werte, maxw) 
   CALL del_params (6, ianz, cpara, lpara, maxw) 
   IF(ianz==1 .OR. ianz==2) THEN
      iianz = 1
      CALL ber_params (ianz, cpara, lpara, wwerte, maxw) 
      IF(ier_num==0) THEN
         cr_spcgrno = NINT(wwerte(1))
         cr_spcgr   = cpara(1)(1:lpara(1))
!        cr_spcgr   = spcgr_name (cr_spcgrno) 
      ELSE
         cr_spcgr = cpara(1)(1:lpara(1))
      ENDIF
      CALL no_error
      CALL del_params (1, ianz, cpara, lpara, maxw) 
      IF(ianz == 1) THEN
         CALL ber_params (ianz, cpara, lpara, wwerte, maxw) 
         IF(ier_num==0) THEN
            spcgr_para = nint (wwerte (1) ) 
            spcgr_ianz = 1
         ELSE
            ier_num = -93
            ier_typ = ER_APPL 
            ier_msg (1) = 'Error reading origin choice indicator'
         ENDIF
      ELSE
         spcgr_para = 1
         spcgr_ianz = 0
      ENDIF
   ENDIF
   DO i = 1, 3 
      cr_a0 (i) = werte (i) 
      cr_win (i) = werte (i + 3) 
   ENDDO 
   IF (cr_a0(1)  <=   0.0 .OR. cr_a0 (2) <=   0.0 .OR. cr_a0 (3) <=   0.0 .OR.    &
       cr_win(1) <=   0.0 .OR. cr_win(2) <=   0.0 .OR. cr_win(3) <=   0.0 .OR.    &
       cr_win(1) >= 180.0 .OR. cr_win(2) >= 180.0 .OR. cr_win(3) >= 180.0      ) THEN                                          
      ier_num = -93 
      ier_typ = ER_APPL 
      ier_msg (1) = 'Error reading unit cell parameters' 
      RETURN                 ! Jump to handle error messages, amd macro conditions
   ENDIF 
   werte(1)=spcgr_para
   iianz = 1
   CALL spcgr_no(iianz,maxw,werte)
ELSE 
   ier_num = - 6 
   ier_typ = ER_COMM 
   RETURN                 ! Jump to handle error messages, amd macro conditions
ENDIF 
IF(cr_syst/=4 .AND. cr_iset/=1) THEN !non-orthorhombic and nonstandard setting
   ier_num = -160
   ier_typ = ER_APPL
   ier_msg(1) = 'Non-standard settings are implemented only'
   ier_msg(2) = 'for orthorhombic space groups '
   RETURN
ENDIF
cr_icc (1) = 1 
cr_icc (2) = 1 
cr_icc (3) = 1 
cr_natoms = 0 
cr_ncatoms = 1 
cr_ncreal  = 1 
cr_nscat = 0 
as_natoms = 0 
!                                                                       
!----reset microdomain status                                      
!                                                                       
CALL do_stack_rese 
!  Flag that no Fourier has been calculated yet
four_last = FOUR_NN
!
chem_purge = .FALSE.    ! No purge was done, period boundary is OK
!
END SUBROUTINE do_readfree
!
!********************************************************************** 
!
SUBROUTINE readcell (strucfile, l_identical, r_identical) 
!-                                                                      
!           This subroutine reads a unit cell.                          
!+                                                                      
use atom_line_mod
USE discus_config_mod 
USE discus_allocate_appl_mod
USE crystal_mod 
USE molecule_mod 
USE prop_para_mod
USE discus_save_mod 
USE spcgr_apply
USE wyckoff_mod
!
use blanks_mod
USE precision_mod
USE lib_errlist_func
USE lib_length
USE str_comp_mod
USE string_convert_mod
USE support_mod
!
IMPLICIT none 
!                                                                       
CHARACTER ( LEN=* ), INTENT(OUT) :: strucfile 
LOGICAL            , INTENT(IN) :: l_identical
REAL               , INTENT(IN) :: r_identical
!
CHARACTER(len=10)          :: befehl 
CHARACTER(LEN=PREC_STRING) :: line, zeile 
INTEGER, PARAMETER                   :: AT_MAXP = 16
INTEGER, PARAMETER                   :: ist     = 7
INTEGER, PARAMETER                   :: MAXW = MAX(AT_MAXP,7)
!
INTEGER                              :: at_ianz
CHARACTER(LEN=8), DIMENSION(AT_MAXP) :: at_param
INTEGER     :: i, j, ibl, lbef 
INTEGER     :: iatom
INTEGER     :: n_read = 0
INTEGER     :: lline 
INTEGER     :: new_nmax
INTEGER     :: new_nscat
INTEGER     :: io_line
INTEGER     :: iimole
INTEGER                          :: n_mole 
INTEGER                          :: n_type 
INTEGER                          :: n_atom 
integer, dimension(3) :: n_cells
LOGICAL          :: need_alloc = .false.
LOGICAL          :: lcontent
LOGICAL, SAVE          :: at_init = .TRUE.
LOGICAL :: lcell, lout 
REAL(KIND=PREC_DP) :: werte (maxw)
REAL :: dw1 , occ1
!                                                                       
LOGICAL :: IS_IOSTAT_END
!                                                                       
cr_natoms = 0 
lcell     = .true. 
lout      = .false. 
lcontent  = .false.
at_param(:) = ' '
at_ianz     = 0
n_read      = 0
CALL test_file(strucfile, new_nmax, new_nscat, n_mole, n_type, &
               n_atom, n_cells, -1 , cr_newtype)
IF (ier_num /= 0) THEN
   CLOSE (ist)
   RETURN
ENDIF
IF( NMAX    < new_nmax .or. &          ! Allocate sufficient atom numbers
    MAXSCAT < new_nscat     ) THEN     ! Allocate sufficient atom types
   new_nmax = MAX(new_nmax ,NMAX)
   new_nscat= MAX(new_nscat,MAXSCAT)
   CALL alloc_crystal(new_nscat, new_nmax)
   IF ( ier_num /= 0) THEN
      CLOSE (IST)
      RETURN
   ENDIF
ENDIF
need_alloc = .false.
IF ( n_mole > MOLE_MAX_MOLE  .or.  &
     n_type > MOLE_MAX_TYPE  .or.  &
     n_atom > MOLE_MAX_ATOM      ) THEN
   n_mole = MAX(n_mole,MOLE_MAX_MOLE,1)
   n_type = MAX(n_type,MOLE_MAX_TYPE,1)
   n_atom = MAX(n_atom,MOLE_MAX_ATOM,1)
   CALL alloc_molecule(1, 1, n_mole, n_type, n_atom)
   IF ( ier_num /= 0) THEN
      CLOSE (IST)
      RETURN
   ENDIF
ENDIF
!
! To ensure that molecules are placed properly, it best to read as
! structure first, shift all molecule origins into the first unit
! cell and to read again
!
IF(n_mole > 0) THEN
   CALL readcell_mole(strucfile, l_identical, r_identical)
   RETURN
ENDIF
!
! No molecules, use old readcell, might be phased out ?????
!
CALL oeffne (ist, strucfile, 'old') 
IF (ier_num /= 0) THEN
   CLOSE (ist)
   RETURN
ENDIF
!
cr_dim(:,1) =  1.e10
cr_dim(:,2) = -1.e10
!     DO i = 1, 3 
!        cr_dim (i, 1) =  1.e10 
!        cr_dim (i, 2) = -1.e10 
!     ENDDO 
!                                                                       
!     --Read header of structure file                                   
!                                                                       
CALL stru_readheader (ist, MAXSCAT, cr_name,      &
         cr_spcgr, cr_set, cr_at_lis, cr_nscat, cr_dw, cr_occ, cr_a0, cr_win, sav_ncell,&
         sav_r_ncell, sav_ncatoms, spcgr_ianz, spcgr_para, AT_MAXP, at_ianz, at_param)
IF (ier_num.ne.0) THEN 
   ier_msg(1) = 'Structure ' // strucfile
   CLOSE (ist)
   RETURN 
ENDIF 
CALL setup_lattice(cr_a0, cr_ar, cr_eps, cr_gten, cr_reps,    &
         cr_rten, cr_win, cr_wrez, cr_v, cr_vr, lout, cr_gmat, cr_fmat, &
         cr_cartesian,                                                  &
         cr_tran_g, cr_tran_gi, cr_tran_f, cr_tran_fi)
!                                                                       
IF (ier_num /= 0) THEN
   CLOSE (ist)
   RETURN
ENDIF
!                                                                       
CALL get_symmetry_matrices 
IF(cr_syst/=4 .AND. cr_iset/=1) THEN !non-orthorhombic and nonstandard setting
   ier_num = -160
   ier_typ = ER_APPL
   ier_msg(1) = 'Non-standard settings are implemented only'
   ier_msg(2) = 'for orthorhombic space groups '
   RETURN
ENDIF
IF( NMAX < spc_n*new_nmax .or.  &      ! Allocate sufficient atom numbers
    MAXSCAT < new_nscat       ) THEN   ! Allocate sufficient scattering types
   new_nmax  = MAX(spc_n*new_nmax + 1, NMAX)
   new_nscat = MAX(new_nscat         , MAXSCAT)
   CALL alloc_crystal(new_nscat, new_nmax)
  IF ( ier_num /= 0) THEN
      CLOSE (IST)
      RETURN
   ENDIF
ENDIF
need_alloc = .false.
IF ( n_mole*spc_n > MOLE_MAX_MOLE ) THEN
   n_mole = n_mole*spc_n
   need_alloc = .true.
ENDIF
IF ( n_type > MOLE_MAX_TYPE ) THEN
   need_alloc = .true.
ENDIF
IF ( n_atom*spc_n > MOLE_MAX_ATOM ) THEN
   n_atom = n_atom*spc_n
   need_alloc = .true.
ENDIF
IF( need_alloc )  THEN         ! Allocate sufficient molecules
   CALL alloc_molecule(1, 1, n_mole, n_type, n_atom)
  IF ( ier_num /= 0) THEN
      CLOSE (IST)
      RETURN
   ENDIF
ENDIF

at_init = .TRUE.

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
   call tab2blank(line, lline)
!23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
empty:   IF (line.ne.' '.and.line (1:1) .ne.'#'.and.line(1:1)/='!' .AND. line.ne.char (13)) THEN
      need_alloc = .false.
      new_nmax   = NMAX
      new_nscat  = MAXSCAT
      IF ( NMAX < cr_natoms + spc_n ) THEN     ! Allocate sufficient atom numbers
         new_nmax  = MAX(NMAX + spc_n + 1, cr_natoms + spc_n+1)
         need_alloc = .true.
      ENDIF
      IF ( MAXSCAT < cr_nscat + 1       ) THEN ! Allocate sufficient atom types
         new_nscat = MAX(MAXSCAT + 5, INT ( MAXSCAT * 1.025 ) )
         need_alloc = .true.
      ENDIF
      IF( need_alloc ) THEN
         CALL alloc_crystal(new_nscat, new_nmax)
         IF ( ier_num /= 0) THEN
            CLOSE (IST)
            RETURN
         ENDIF
         ier_num = -49 
         ier_typ = ER_APPL 
      ENDIF
!
      lbef = 10 
      befehl = ' ' 
      ibl = index (line (1:lline) , ' ') 
      IF (ibl.eq.0) THEN 
         ibl = lline+1 
      ENDIF 
      lbef = min (ibl - 1, lbef) 
      befehl = line (1:lbef) 
typus:IF(str_comp (befehl, 'molecule', 4, lbef, 8) .or.       &
         str_comp (befehl, 'domain',   4, lbef, 6) .or.       &
         str_comp (befehl, 'object',   4, lbef, 6)     ) THEN
!                                                                       
!     ----------Start/End of a molecule                                 
!                                                                       
         CALL no_error 
         IF(ibl.le.lline) THEN 
            i = lline-ibl 
            zeile = line (ibl + 1:lline) 
         ELSE 
            zeile = ' ' 
            i = 0 
         ENDIF 
         CALL struc_mole_header (zeile, i, .true., lcontent) 
         IF (ier_num.ne.0) THEN
            CLOSE(IST)
            RETURN 
         ENDIF
      ELSE  typus
         DO j = 1, MAXW 
            werte (j) = 0.0 
         ENDDO 
         werte (5) = 1.0 
         cr_surf(:,cr_natoms+1) = 0
         cr_magn(:,cr_natoms+1) = 0.0
         as_natoms = as_natoms + 1 
         call get_atom_werte(as_natoms, MAXW, werte)
!                  CALL read_atom_line (line, ibl, lline, as_natoms, MAXW, werte, &
!                                       AT_MAXP, at_ianz, at_param, at_init)                                          
         n_read = n_read + 1
         IF (ier_num.ne.0.and.ier_num.ne. -49) THEN 
            GOTO 999 
         ENDIF 
         cr_natoms = cr_natoms + 1 
         i = cr_natoms 
         IF (.not.mole_l_on .AND. NINT(werte(6))==0 ) THEN 
!                                                                       
!     ------------Transform atom into first unit cell,                  
!                 if it is not inside a molecule                        
!                                                                       
            CALL firstcell (werte, maxw) 
         ENDIF 
         IF(NINT(werte(6))>0 .AND. NINT(werte(7))>0) THEN
            iimole = NINT(werte(6)) !+ n_mole_old
            CALL mole_insert_explicit(cr_natoms, iimole        , NINT(werte(7))) 
            cr_prop(cr_natoms) = IBSET(cr_prop(cr_natoms),PROP_MOLECULE)
            cr_mole(cr_natoms) = NINT(werte(6)) !+ n_mole_old
!                  ELSEIF(mole_l_on) THEN
!                     CALL mole_insert_current (cr_natoms, mole_num_curr) 
         ENDIF 
!                                                                       
         DO j = 1, 3 
            cr_pos (j, i) = werte (j) 
         ENDDO 
         cr_surf(0:3,i) = NINT(werte(9:12))
         cr_magn(0:3,i) = 0.0 ! (werte(9:12)) MAGNETIC_WORK
         dw1 = werte (4) 
         occ1 = werte(8)                       ! WORK OCC
         IF(mole_l_on) THEN
!           cr_mole (i) = mole_num_mole
         ELSE
            cr_mole(i) = NINT(werte(6))
         ENDIF
         cr_prop (i) = NINT(werte(5) ) 
         cr_surf(:,i) = 0                      ! Currently no save nor read for surface
         cr_magn(:,i) = werte(13:16)           ! Read for MAGNETIC_WORK
         IF(MAXVAL(ABS(werte(13:16)))> 0.0) cr_magnetic = .TRUE.  ! Crystal has magnetic atoms
!                                                                       
         if_blk: IF(line(1:ibl) .ne.'    ') THEN 
            ibl = ibl - 1 
            CALL do_cap (line (1:ibl) ) 
!                                                                       
!------ ----------- New option determines whether all atoms in          
!------ ----------- asymmetric unit are considered different atom       
!------ ----------- types ..                                            
!                                                                       
            IF (.not.cr_newtype) THEN          ! lcell == .true.
               DO j = 0, cr_nscat 
                  IF (line (1:ibl)  == cr_at_lis (j)              &
                        .and.dw1 == cr_dw (j)                     &
                        .AND. occ1==cr_occ(j) ) THEN                    
                     cr_iscat (i) = j 
                     CALL symmetry 
                     IF (ier_num.ne.0) THEN 
                        CLOSE (IST)
                        RETURN 
                     ENDIF 
                     GOTO 22 
!                       ELSEIF(line(1:ibl)=='VOID' .AND. mole_l_on) THEN
!                          cr_iscat (i) = 0 
!                          CALL symmetry 
!                          IF (ier_num.ne.0) THEN 
!                             CLOSE (IST)
!                             RETURN 
!                          ENDIF 
!                          GOTO 22 
                  ENDIF 
               ENDDO 
            ENDIF 
!                                                                       
!------ ----------- end new code                                        
!                                                                       
            IF (cr_nscat.lt.maxscat) THEN 
!ATOM_LINE              as_natoms = as_natoms + 1 
               cr_nscat = cr_nscat + 1 
               cr_iscat (i) = cr_nscat 
               cr_at_lis (cr_nscat) = line (1:ibl) 
               cr_dw (cr_nscat) = dw1 
               cr_occ(cr_nscat) = occ1                    ! WORK OCC
!                                                                       
!              as_at_lis (cr_nscat) = cr_at_lis (cr_nscat) 
!              as_iscat (as_natoms) = cr_iscat (i) 
!              as_dw (as_natoms) = cr_dw (cr_nscat) 
!              as_occ(as_natoms) = cr_occ(cr_nscat) 
!                       ENDIF
!              DO j = 1, 3 
!                 as_pos (j, as_natoms) = cr_pos (j, i) 
!              ENDDO 
!              as_mole (as_natoms) = cr_mole (i) 
!              as_prop (as_natoms) = cr_prop (i) 
               CALL symmetry 
               IF (ier_num.ne.0) THEN 
               CLOSE(IST)
               RETURN 
               ENDIF 
            ELSE 
               ier_num = -26 
               ier_typ = ER_APPL 
               GOTO 2 
            ENDIF 
   22       CONTINUE 
         ENDIF  if_blk
      ENDIF  typus
   ENDIF empty
ENDDO main 
!
CALL test_identical (l_identical, r_identical) ! Test if atoms are too close
!                                                                       
2    CONTINUE 
IF (ier_num.eq. -49) THEN 
   CALL no_error 
!                                                                       
!       move first unit cell into lower left corner of crystal          
!                                                                       
   DO i = 1, cr_natoms 
      DO j = 1, 3 
         cr_pos (j, i) = cr_pos (j, i) - int ( (cr_icc (j) ) / 2) 
!             cr_pos(j,i)=cr_pos(j,i) - int((cr_icc(j)-0.1)/2)          
         cr_dim (j, 1) = min (cr_dim (j, 1), cr_pos (j, i) ) 
         cr_dim (j, 2) = max (cr_dim (j, 2), cr_pos (j, i) ) 
      ENDDO 
   ENDDO 
ELSE
   CLOSE(IST)
   WRITE (ier_msg (1), 3000) n_read
   RETURN
ENDIF 
!     INTEGER     :: iatom
!
!     If a molecule containted a "molecule atoms" instruction, we need to
!     set the molecule flag
!
IF(lcontent) THEN 
   DO i = 1, mole_num_mole
      DO j = 1, mole_len (i)
         iatom          = mole_cont (mole_off (i) + j)
         cr_prop(iatom) = ibset(cr_prop(iatom),PROP_MOLECULE)
         cr_mole(iatom) = i
      ENDDO
   ENDDO
ENDIF 
!                                                                       
!     ENDIF 
!                                                                       
  999 CONTINUE 
CLOSE (ist) 
IF (ier_num.eq. - 49) THEN 
   WRITE (ier_msg (1), 3000) as_natoms + 1 
 3000 FORMAT      ('At atom number = ',i8) 
ENDIF 
!
CALL test_mole_gap
!                                                                       
 2000 FORMAT    (a) 
!
END SUBROUTINE readcell                       
!
!********************************************************************** 
!
      SUBROUTINE struc_mole_header (zeile, lp, lcell, lcontent) 
!-                                                                      
!     interprets the 'molecule' lines of a structure file               
!+                                                                      
                                                                        
      USE discus_allocate_appl_mod
      USE discus_config_mod 
      USE crystal_mod 
      USE molecule_mod 
      USE spcgr_apply
      USE ber_params_mod
      USE get_params_mod
USE precision_mod
USE str_comp_mod
      IMPLICIT none 
!                                                                       
      CHARACTER(LEN=* ), INTENT(IN)    :: zeile 
      INTEGER          , INTENT(INOUT) :: lp
      LOGICAL          , INTENT(IN)    :: lcell 
      LOGICAL          , INTENT(OUT)   :: lcontent 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 21) 
!                                                                       
      CHARACTER(LEN=MAX(PREC_STRING,LEN(zeile))) :: cpara (maxw) 
      INTEGER j, ianz 
      INTEGER lpara (maxw)
      REAL(KIND=PREC_DP) :: werte (maxw) 
      INTEGER          :: n_gene
      INTEGER          :: n_symm
      INTEGER          :: n_mole
      INTEGER          :: n_type
      INTEGER          :: n_atom
      LOGICAL          :: need_alloc = .false.
!                                                                       
!
!
CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
!
IF (ier_num.eq.0) THEN 
   IF (ianz.eq.0) THEN 
!                                                                       
!     --No parameters, start a new Molekule                             
!                                                                       
      need_alloc = .false.
      n_gene = MAX( 1, MOLE_MAX_GENE)
      n_symm = MAX( 1, MOLE_MAX_SYMM)
      n_mole =         MOLE_MAX_MOLE
      n_type =         MOLE_MAX_TYPE
      n_atom =         MOLE_MAX_ATOM
      IF (mole_num_mole >= MOLE_MAX_MOLE ) THEN
         n_mole = mole_num_mole + 20
         need_alloc = .true.
      ENDIF
      IF (mole_num_type >= MOLE_MAX_TYPE ) THEN
         n_type = mole_num_type + 10
         need_alloc = .true.
      ENDIF
      IF ( need_alloc ) THEN
         call alloc_molecule(n_gene, n_symm, n_mole, n_type, n_atom)
      ENDIF
      IF (mole_num_mole.lt.MOLE_MAX_MOLE) THEN 
         IF (mole_num_type.lt.MOLE_MAX_TYPE) THEN 
            mole_l_on = .true. 
            mole_l_first = .true. 
            mole_num_atom = mole_off (mole_num_mole) + &
                            mole_len (mole_num_mole)                                        
            mole_num_mole = mole_num_mole+1 
            mole_num_curr = mole_num_mole 
            mole_num_type = mole_num_type+1 
            mole_off (mole_num_mole) = mole_num_atom 
            mole_type(mole_num_mole) = mole_num_type 
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
      IF (str_comp (cpara (1) , 'end', 3, lpara (1) , 3) ) THEN 
!                                                                       
!     ----Turn off molecule                                             
!                                                                       
         IF(mole_len(mole_num_mole)>0) THEN    !Ignore empty molecules
            IF (lcell) call mole_firstcell 
         ENDIF
           mole_l_on = .false. 
!                                                                       
      ELSEIF(str_comp(cpara(1), 'character', 3, lpara(1), 9)) THEN
!                                                                       
!     ------Define whether this is a molecule or an object              
!                                                                       
         IF(str_comp(cpara(2), 'atoms', 2, lpara(2), 5) ) THEN
            mole_char (mole_num_mole) = MOLE_ATOM 
         ELSEIF (str_comp (cpara (2) , 'cube', 2, lpara (2) , 4) ) THEN
            mole_char (mole_num_mole) = MOLE_CUBE 
         ELSEIF (str_comp (cpara (2) , 'cylinder', 2, lpara (2) , 8) ) THEN
            mole_char (mole_num_mole) = MOLE_CYLINDER 
         ELSEIF (str_comp (cpara (2) , 'sphere', 2, lpara (2) , 6)) THEN
            mole_char (mole_num_mole) = MOLE_SPHERE 
         ELSEIF (str_comp (cpara (2) , 'edge', 2, lpara (2) , 4) ) THEN
            mole_char (mole_num_mole) = MOLE_EDGE 
         ELSEIF (str_comp (cpara (2) , 'domain_cube', 9, lpara (2), 11) ) THEN
            mole_char (mole_num_mole) = MOLE_DOM_CUBE 
         ELSEIF (str_comp (cpara (2) , 'domain_cylinder', 9, lpara (2), 15) ) THEN
            mole_char (mole_num_mole) = MOLE_DOM_CYLINDER 
         ELSEIF (str_comp (cpara (2) , 'domain_sphere', 9, lpara(2), 13) ) THEN
            mole_char (mole_num_mole) = MOLE_DOM_SPHERE 
         ELSEIF (str_comp (cpara (2) , 'domain_fuzzy', 9, lpara(2), 12) ) THEN
            mole_char (mole_num_mole) = MOLE_DOM_FUZZY 
         ELSE 
            ier_num = - 82 
            ier_typ = ER_APPL 
         ENDIF 
!                                                                       
      ELSEIF (str_comp (cpara (1) , 'file', 3, lpara (1) , 4) ) THEN
         mole_file (mole_num_mole) = cpara (2) (1:lpara(2))
!                                                                       
      ELSEIF (str_comp (cpara (1) , 'density', 3, lpara (1) , 6) ) THEN
!                                                                       
!     ------Define the scattering density of an object                  
!                                                                       
         cpara (1) = '0' 
         lpara (1) = 1 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.eq.0) THEN 
            IF (ianz.eq.2) THEN 
               mole_dens (mole_num_mole) = werte (2) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
         ENDIF 
!                                                                       
      ELSEIF(str_comp(cpara(1),'biso',3,lpara(1),3)) THEN                                                        
!                                                                       
!     ------Define the isotropic molecular B-Value
!                                                                       
         cpara (1) = '0' 
         lpara (1) = 1 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.eq.0) THEN 
            IF (ianz.eq.2) THEN 
               mole_biso(mole_type(mole_num_mole)) = werte (2) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
         ENDIF 
!                                                                       
      ELSEIF(str_comp(cpara(1),'clin',3,lpara(1),3)) THEN                                                        
!                                                                       
!     ------Define the isotropic molecular B-Value
!                                                                       
         cpara (1) = '0' 
         lpara (1) = 1 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.eq.0) THEN 
            IF (ianz.eq.2) THEN 
               mole_clin(mole_type(mole_num_mole)) = werte (2) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
         ENDIF 
!                                                                       
      ELSEIF(str_comp(cpara(1),'cqua',3,lpara(1),3)) THEN                                                        
!                                                                       
!     ------Define the isotropic molecular B-Value
!                                                                       
         cpara (1) = '0' 
         lpara (1) = 1 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.eq.0) THEN 
            IF (ianz.eq.2) THEN 
               mole_cqua(mole_type(mole_num_mole)) = werte (2) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
         ENDIF 
!                                                                       
      ELSEIF (str_comp (cpara (1) , 'fuzzy', 3, lpara (1) , 5) ) THEN
!                                                                       
!     ------Define the minimum distance between atoms in a domain       
!             and the host                                              
!                                                                       
         cpara (1) = '0' 
         lpara (1) = 1 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.eq.0) THEN 
            IF (ianz.eq.2) THEN 
               mole_fuzzy (mole_num_mole) = werte (2) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
         ENDIF 
!                                                                       
      ELSEIF (str_comp (cpara (1) , 'generator', 3, lpara (1) , 9)) THEN
!                                                                       
!     ------Define which generators create atoms within the             
!           same molecule                                               
!           Obsolete statement, is done automatically!!! RBN            
!                                                                       
         ier_num = + 2 
         ier_typ = ER_APPL 
!                                                                       
      ELSEIF (str_comp (cpara (1) , 'symmetry', 4, lpara (1) , 8) ) THEN
!                                                                       
!     ------Define which symmetries  create atoms within the            
!           same molecule                                               
!           Obsolete statement, is done automatically!!! RBN            
!                                                                       
         ier_num = + 2 
         ier_typ = ER_APPL 
!                                                                       
      ELSEIF (str_comp (cpara (1) , 'type', 3, lpara (1) , 4) ) THEN
!                                                                       
!     ------Define the molecule type, if less than current highest      
!     ------type number diminuish type number by one                    
!                                                                       
         cpara (1) = '0' 
         lpara (1) = 1 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.eq.0) THEN 
            IF (ianz.eq.2) THEN 
               IF (NINT (werte (2) ) .lt.mole_num_type) THEN 
                  mole_num_type = mole_num_type-1 
                  mole_type (mole_num_mole) = NINT (werte (2) ) 
               ELSE
                  mole_type (mole_num_mole) = NINT (werte (2) )
                  mole_num_type = MAX(mole_num_type, mole_type (mole_num_mole))
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
         ENDIF 
!                                                                       
      ELSEIF (str_comp (cpara (1) , 'content', 4, lpara (1) , 7) ) THEN
!                                                                       
!     ------start reading a molecule content                            
!                                                                       
         cpara (1) = '0' 
         lpara (1) = 1 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.eq.0) THEN 
            IF (ianz.eq.2.or.ianz.eq.3) THEN 
               IF (mole_num_mole.lt.MOLE_MAX_MOLE) THEN 
                  IF (werte (2) .lt.MOLE_MAX_TYPE) THEN 
                     IF (mole_l_on) THEN 
                        mole_type(mole_num_mole) = int(werte(2))
                        mole_num_type = max (mole_num_type-1,     &
                                                 int(werte (2) ) )                            
                    ELSE 
                        mole_num_atom = mole_off (mole_num_mole)  &
                            + mole_len (mole_num_mole)                
                        mole_num_mole = mole_num_mole+1 
                        mole_num_curr = mole_num_mole 
                        mole_type (mole_num_mole) = int (werte (2))
                        mole_off (mole_num_mole) = mole_num_atom 
                        mole_len (mole_num_mole) = 0 
                        mole_num_type = max (mole_num_type,       &
                                                 int(werte (2) ) )                             
                        mole_gene_n = 0 
                        mole_symm_n = 0 
                        mole_l_on = .true. 
                    ENDIF 
                  ELSE 
                     ier_num = - 64 
                     ier_typ = ER_APPL 
                     ier_msg (1) = 'First characters of wrong line' 
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
      ELSEIF (str_comp (cpara (1) , 'atoms', 4, lpara (1) , 5) ) THEN
!                                                                       
!     ------read a molecule content                                     
!                                                                       
         IF (mole_l_on) THEN 
            cpara (1) = '0' 
            lpara (1) = 1 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            IF (ier_num.eq.0) THEN 
               DO j = 2, ianz 
                  mole_len (mole_num_mole) = mole_len (mole_num_mole)+1
                  mole_cont (mole_off (mole_num_mole) + &
                             mole_len ( mole_num_mole) ) = int (werte (j) )                
               ENDDO 
               lcontent = .true.
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
      cr_set,                                                           &
      cr_a0, cr_win, cr_natoms, cr_nscat, cr_dw, cr_occ, cr_at_lis, cr_pos,     &
      cr_mole, cr_surf, cr_magn,                                        &
      cr_iscat, cr_prop, cr_dim, cr_magnetic,                           &
      as_natoms, as_at_lis, as_dw, as_pos,   &
      as_iscat, as_prop, sav_ncell, sav_r_ncell, sav_ncatoms,           &
      spcgr_ianz, spcgr_para)                                           
!-                                                                      
!           this subroutine reads an old structur.                      
!+                                                                      
!
use precision_mod
USE support_mod
!
IMPLICIT none 
!                                                                       
INTEGER                                  , intent(in)    :: NMAX 
INTEGER                                  , intent(inout) :: MAXSCAT 
CHARACTER(len=*)                         , intent(in)    :: strucfile 
CHARACTER(len=80)                        , intent(out)   :: cr_name 
CHARACTER(len=16)                        , intent(out)   :: cr_spcgr 
CHARACTER(LEN=3)                         , INTENT(INOUT) :: cr_set
REAL(kind=PREC_DP), dimension(3)         , intent(out)   :: cr_a0 (3) 
REAL(kind=PREC_DP), dimension(3)         , intent(out)   :: cr_win (3) 
INTEGER                                  , INTENT(INOUT) :: cr_natoms
INTEGER                                  , intent(inout) :: cr_nscat 
REAL(kind=PREC_DP), dimension(0:MAXSCAT) , intent(out  ) :: cr_dw ! (0:MAXSCAT) 
REAL(kind=PREC_DP), dimension(0:MAXSCAT) , intent(out  ) :: cr_occ ! (0:MAXSCAT) 
CHARACTER(len=4)  , dimension(0:MAXSCAT) , intent(out  ) :: cr_at_lis ! (0:MAXSCAT) 
REAL(kind=PREC_DP), DIMENSION(1:3,1:NMAX), INTENT(out  ) :: cr_pos
INTEGER           , DIMENSION(1:NMAX),     INTENT(out  ) :: cr_mole
INTEGER           , DIMENSION(0:3,1:NMAX), INTENT(out  ) :: cr_surf
REAL(kind=PREC_DP), DIMENSION(0:3,1:NMAX), INTENT(out  ) :: cr_magn
INTEGER           , DIMENSION(1:NMAX),     INTENT(out  ) :: cr_iscat
INTEGER           , DIMENSION(1:NMAX),     INTENT(out  ) :: cr_prop
REAL(kind=PREC_DP), dimension(3,2)       , intent(out  ) :: cr_dim ! (3, 2) 
LOGICAL                                  , INTENT(INOUT) :: cr_magnetic
INTEGER                                  , intent(inout) :: as_natoms 
CHARACTER(len=4)  , dimension(0:MAXSCAT) , intent(inout) :: as_at_lis ! (0:MAXSCAT) 
REAL(kind=PREC_DP), dimension(0:MAXSCAT) , intent(out)   :: as_dw ! (0:MAXSCAT) 
REAL(kind=PREC_DP), dimension(3, MAXSCAT), intent(out)   :: as_pos ! (3, MAXSCAT) 
INTEGER           , dimension(MAXSCAT)   , intent(out)   :: as_iscat ! (MAXSCAT) 
INTEGER           , dimension(MAXSCAT)   , intent(out)   :: as_prop  !(MAXSCAT) 
INTEGER                                  , intent(inout) :: sav_ncell (3) 
LOGICAL                                  , intent(inout) :: sav_r_ncell 
INTEGER                                  , intent(inout) :: sav_ncatoms 
INTEGER                                  , intent(out)   :: spcgr_ianz 
INTEGER                                  , intent(out)   :: spcgr_para 
!
!                                                                       
INTEGER, PARAMETER                   :: ist = 7 
INTEGER, PARAMETER                   :: AT_MAXP = 16
!
INTEGER :: i 
INTEGER                              :: at_ianz
LOGICAL :: lcell 
CHARACTER(LEN=8), DIMENSION(AT_MAXP) :: at_param
!                                                                       
REAL(kind=PREC_DP), dimension(0:MAXSCAT) :: as_occ(0:MAXSCAT) 
!                                                                       
!                                                                       
      cr_natoms = 0 
      lcell = .false. 
      at_param(:) = ' '
      at_ianz     = 0
      CALL oeffne (ist, strucfile, 'old') 
      IF (ier_num.eq.0) THEN 
         DO i = 1, 3 
         cr_dim (i, 1) = 1.D10 
         cr_dim (i, 2) = - 1.D10 
         ENDDO 
!                                                                       
!     --Read header of structure file                                   
!                                                                       
         CALL stru_readheader (ist, MAXSCAT, cr_name,      &
         cr_spcgr, cr_set, cr_at_lis, cr_nscat, cr_dw, cr_occ, cr_a0, cr_win, sav_ncell,&
         sav_r_ncell, sav_ncatoms, spcgr_ianz, spcgr_para, AT_MAXP, at_ianz, at_param)
!                                                                       
         IF (ier_num.eq.0) THEN 
!                                                                       
            CALL struc_read_atoms (NMAX, MAXSCAT, cr_natoms, cr_nscat,  &
            cr_dw, cr_occ, cr_at_lis, cr_pos, cr_iscat, cr_mole, cr_surf, &
            cr_magn, cr_prop, cr_dim, cr_magnetic, &
            as_natoms, as_at_lis, as_dw, as_occ, as_pos, as_iscat, as_prop, &
            AT_MAXP, at_ianz, at_param)     
         else
            ier_msg(1) = 'Structure ' // strucfile
         ENDIF 
      ENDIF 
!                                                                       
      CLOSE (ist) 
      IF (ier_num.eq. -49) THEN 
         WRITE (ier_msg (1), 3000) cr_natoms + 1 
         ier_msg(2) = 'Structure ' // strucfile
 3000 FORMAT      ('At atom number = ',i8) 
      ENDIF 
      END SUBROUTINE readstru                       
!********************************************************************** 
SUBROUTINE stru_readheader (ist, HD_MAXSCAT, cr_name,   &
      cr_spcgr, cr_set, cr_at_lis, cr_nscat, cr_dw, cr_occ, cr_a0, cr_win, sav_ncell,   &
      sav_r_ncell, sav_ncatoms, spcgr_ianz, spcgr_para, AT_MAXP, at_ianz, at_param)                 
!-                                                                      
!     This subroutine reads the header of a structure file              
!+                                                                      
USE gen_add_mod 
USE sym_add_mod 
!
use blanks_mod
USE ber_params_mod
USE get_params_mod
USE lib_errlist_func
USE lib_length
USE string_convert_mod
USE precision_mod
USE take_param_mod
USE str_comp_mod
!
IMPLICIT none 
!                                                                       
INTEGER                                  , INTENT(IN)  :: ist
INTEGER                                  , INTENT(IN)  :: HD_MAXSCAT 
CHARACTER(LEN=80)                        , INTENT(OUT) :: cr_name 
CHARACTER(LEN=16)                        , INTENT(OUT) :: cr_spcgr 
CHARACTER(LEN= 3)                        , INTENT(OUT) :: cr_set 
CHARACTER(LEN=4), DIMENSION(0:HD_MAXSCAT), INTENT(OUT) :: cr_at_lis ! (0:HD_MAXSCAT) 
INTEGER                                  , INTENT(OUT) :: cr_nscat 
REAL(kind=PREC_DP), DIMENSION(0:HD_MAXSCAT)         , INTENT(OUT) :: cr_dw     ! (0:HD_MAXSCAT) 
REAL(kind=PREC_DP), DIMENSION(0:HD_MAXSCAT)         , INTENT(OUT) :: cr_occ    ! (0:HD_MAXSCAT) 
REAL(kind=PREC_DP), DIMENSION(3)                    , INTENT(OUT) :: cr_a0     ! (3) 
REAL(kind=PREC_DP), DIMENSION(3)                    , INTENT(OUT) :: cr_win    ! (3) 
INTEGER, DIMENSION(3)                    , INTENT(OUT) :: sav_ncell ! (3) 
INTEGER                                  , INTENT(OUT) :: sav_ncatoms 
LOGICAL                                  , INTENT(OUT) :: sav_r_ncell 
!                                                                       
INTEGER                                  , INTENT(OUT) :: spcgr_ianz 
INTEGER                                  , INTENT(OUT) :: spcgr_para 
INTEGER                                  , INTENT(IN)  :: AT_MAXP
INTEGER                                  , INTENT(OUT) :: at_ianz
CHARACTER(LEN=8), DIMENSION(AT_MAXP)     , INTENT(OUT) :: at_param
!                                                                       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 13) 
!                                                                       
      CHARACTER(LEN=PREC_STRING) :: line, cpara (maxw) 
      CHARACTER(LEN=PREC_STRING) :: zeile 
      CHARACTER(LEN=6) :: befehl 
      INTEGER i, ll, j, islash 
      INTEGER ianz 
!DBG      integer             spcgr_ianz                                
      INTEGER lpara (maxw), lp 
      INTEGER lbef, indxb 
      INTEGER xx_nscat, xx_nadp , xx_nocc
      LOGICAL lend 
      LOGICAL :: lcontent
!DBG      real            spcgr_para                                    
      REAL(KIND=PREC_DP) :: werte (maxw) 
!
INTEGER, PARAMETER :: NOPTIONAL = 1
INTEGER, PARAMETER :: O_SETTING = 1
CHARACTER(LEN=   7)       , DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=PREC_STRING), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 0 ! Number of values to calculate 
!                                                                       
!                                                                       
!
DATA oname  / 'setting' /
DATA loname /  7        /
!
opara  =  (/ 'abc' /)   ! Always provide fresh default values
lopara =  (/  3   /)
owerte =  (/  0.0 /)
!
cr_occ(:) = 1.0D0  !! WORK OCC
      xx_nscat = 0 
      xx_nadp = 0 
      xx_nocc = 0 
gen_add_n = 0
sym_add_n = 0
!                                                                       
      ier_num = -46 
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
      IF (str_comp (befehl, 'title', 1, lbef, 5) ) THEN 
!                                                                       
!     Read new header                                                   
!                                                                       
!                                                                       
!     --remove title string from cr_name                                
!                                                                       
         ll = len (cr_name) 
         ll = len_str (cr_name) 
         line = ' ' 
         IF (0.lt.indxb.and.indxb + 1.le.ll) THEN 
            line (1:ll - indxb) = cr_name (indxb + 1:ll) 
         ELSE 
            line = ' ' 
         ENDIF 
         cr_name = line (1:len(cr_name))
!                                                                       
         CALL no_error 
         DO while (.not.str_comp (befehl, 'atoms', 1, lbef, 5) ) 
         READ (ist, 2000, end = 999, err = 999) line 
         ll = 200 
         ll = len_str (line) 
         call tab2blank(line, ll)
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
         IF (indxb + 1.le.ll) THEN 
            zeile = line (indxb + 1:ll) 
            lp = ll - indxb 
         ENDIF 
!                                                                       
!     ----Commentary                                                    
!                                                                       
                                                                        
         IF (line.eq.' '.or.line (1:1) .eq.'#'.or. line(1:1) == '!' .OR. &
             line.eq.char (13) ) THEN
            CONTINUE 
!                                                                       
!     ----Space group symbol                                            
!                                                                       
         ELSEIF (str_comp (befehl, 'spcgr', 1, lbef, 5) ) THEN 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ianz.lt.1) THEN 
               ier_num = - 100 
               ier_typ = ER_APPL 
               RETURN
            ENDIF 
!           Optionally get setting
            CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                              oname, loname, opara, lopara, lpresent, owerte)
            cr_set = opara(O_SETTING)(1:MIN(3,lopara(O_SETTING)))
            islash = index (cpara (1) (1:lpara (1) ) , 'S') 
            DO while (islash.ne.0) 
            cpara (1) (islash:islash) = '/' 
            islash = index (cpara (1) (1:lpara (1) ) , 'S') 
            ENDDO 
            cr_spcgr = cpara (1) (1:lpara(1))
            spcgr_ianz = ianz - 1 
            ianz = ianz - 1 
            spcgr_para = 1 
            IF (ianz.eq.1) THEN 
               cpara (1) = cpara (2) 
               lpara (1) = lpara (2) 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               spcgr_para = nint (werte (1) ) 
            ENDIF 
            IF (ier_num.ne.0) THEN 
               ier_num = - 47 
               ier_typ = ER_APPL 
               RETURN
            ENDIF 
!                                                                       
!     ----Cell constants                                                
!                                                                       
         ELSEIF (str_comp (befehl, 'cell', 1, lbef, 4) ) THEN 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ier_num.eq.0) THEN 
               IF (ianz.eq.6) THEN 
!     --------New style, kommata included                               
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  IF (ier_num.eq.0) THEN 
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
            THEN                                                        
               ier_num = - 93 
               ier_typ = ER_APPL 
               ier_msg (1) = 'Error reading unit cell parameters' 
               RETURN
            ENDIF 
!                                                                       
!     ----Additional symmetry generators 'generator'                    
!                                                                       
         ELSEIF (str_comp (befehl, 'gene', 1, lbef, 4) ) THEN 
            IF (gen_add_n.lt.GEN_ADD_MAX) THEN 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ier_num.eq.0) THEN 
                  IF (ianz.eq.12.or.ianz.eq.13) THEN 
                     CALL ber_params (ianz, cpara, lpara, werte,     &
                     maxw)                                           
                     IF (ier_num.eq.0) THEN 
                        gen_add_n = gen_add_n + 1 
                        DO j = 1, 4 
                        gen_add (1, j, gen_add_n) = werte (j) 
                        gen_add (2, j, gen_add_n) = werte (j + 4) 
                        gen_add (3, j, gen_add_n) = werte (j + 8) 
                        ENDDO 
                        IF (ianz.eq.13) THEN 
                           gen_add_power (gen_add_n) = nint (werte ( &
                           13) )                                     
                        ELSE 
                           gen_add_power (gen_add_n) = 1 
                        ENDIF 
                     ENDIF 
                  ELSEIF (ianz.eq.1) THEN 
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
                     IF (lend) THEN 
                        ier_num = - 92 
                        ier_typ = ER_APPL 
                        RETURN
                     ENDIF 
                  ELSE 
                     ier_num = - 92 
                     ier_typ = ER_APPL 
                     RETURN
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
         ELSEIF (str_comp (befehl, 'scat', 2, lbef, 4) ) THEN 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ier_num.eq.0) THEN 
               IF (xx_nscat + ianz.le.HD_MAXSCAT) THEN 
                  DO i = 1, ianz 
                  CALL do_cap (cpara (i) (1:lpara (i) ) ) 
                  cr_at_lis (xx_nscat + i) = cpara (i) (1:lpara(i))
                  ENDDO 
                  xx_nscat = xx_nscat + ianz 
                  cr_nscat = max (cr_nscat, xx_nscat) 
               ELSE 
                  ier_num = -26 
                  ier_typ = ER_APPL 
                  RETURN
               ENDIF 
            ELSE 
               ier_num = -111
               ier_typ = ER_APPL 
               CLOSE(99)
               RETURN
            ENDIF 
!                                                                       
!     ----Displacement parameters to setup specific sequence of         
!                                    scattering curves 'adp'            
!                                                                       
         ELSEIF (str_comp (befehl, 'adp', 2, lbef, 3) ) THEN 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ier_num.eq.0) THEN 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.eq.0) THEN 
                  IF (xx_nadp + ianz.le.HD_MAXSCAT) THEN 
                     DO i = 1, ianz 
                     cr_dw (xx_nadp + i) = werte (i) 
                     ENDDO 
                     xx_nadp = xx_nadp + ianz 
                     cr_nscat = max (cr_nscat, xx_nadp) 
                  ELSE 
                     ier_num = - 26 
                     ier_typ = ER_APPL 
                     RETURN
                  ENDIF 
               ENDIF 
            ELSE 
               ier_num = -112
               ier_typ = ER_COMM 
               RETURN
            ENDIF 
!                                                                       
!     ----Occupancy parameters to setup specific sequence of         
!                                    scattering curves 'occ'            
!                                                                       
         ELSEIF (str_comp (befehl, 'occ', 2, lbef, 3) ) THEN 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ier_num.eq.0) THEN 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.eq.0) THEN 
                  IF (xx_nocc + ianz.le.HD_MAXSCAT) THEN 
                     DO i = 1, ianz 
                     cr_occ(xx_nocc + i) = werte (i) 
                     ENDDO 
                     xx_nocc = xx_nocc + ianz 
                     cr_nscat = max (cr_nscat, xx_nocc) 
                  ELSE 
                     ier_num = - 26 
                     ier_typ = ER_APPL 
                     RETURN
                  ENDIF 
               ENDIF 
            ELSE 
               ier_num = -112
               ier_typ = ER_COMM 
               RETURN
            ENDIF 
!                                                                       
!     ----Crystal dimensions and number of atoms per unit cell 'ncell'  
!                                                                       
         ELSEIF (str_comp (befehl, 'ncell', 1, lbef, 5) ) THEN 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ier_num.eq.0) THEN 
               IF (ianz.eq.4 .or. ianz==5 .OR. ianz==9) THEN    ! allow for number of atoms on ncell command
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  IF (ier_num.eq.0) THEN 
                     DO j = 1, 3 
                     sav_ncell (j) = int( werte (j) )
                     ENDDO 
                     sav_ncatoms = int( werte (4) )
                     sav_r_ncell = .true. 
                  ENDIF 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
                  RETURN
               ENDIF 
            ENDIF 
!                                                                       
!     ----Additional symmetry operations 'symmetry'                     
!                                                                       
         ELSEIF (str_comp (befehl, 'symm', 2, lbef, 4) ) THEN 
            IF (sym_add_n.lt.SYM_ADD_MAX) THEN 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ier_num.eq.0) THEN 
                  IF (ianz.eq.12.or.ianz.eq.13) THEN 
                     CALL ber_params (ianz, cpara, lpara, werte,     &
                     maxw)                                           
                     IF (ier_num.eq.0) THEN 
                        sym_add_n = sym_add_n + 1 
                        DO j = 1, 4 
                        sym_add (1, j, sym_add_n) = werte (j) 
                        sym_add (2, j, sym_add_n) = werte (j + 4) 
                        sym_add (3, j, sym_add_n) = werte (j + 8) 
                        ENDDO 
                        IF (ianz.eq.13) THEN 
                           sym_add_power (sym_add_n) = nint (werte ( &
                           13) )                                     
                        ELSE 
                           sym_add_power (sym_add_n) = 1 
                        ENDIF 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                     RETURN
                  ENDIF 
               ENDIF 
            ELSE 
               ier_num = - 62 
               ier_typ = ER_APPL 
               RETURN
            ENDIF 
!                                                                       
         ELSEIF (str_comp (befehl, 'molecule', 2, lbef, 8) .OR. & 
                 str_comp (befehl, 'domain',   4, lbef, 6) .OR. &
                 str_comp (befehl, 'object',   4, lbef, 6)     ) THEN
!                                                                       
!     ------Start/End of a molecule                                     
!                                                                       
            CALL no_error 
            IF (indxb.le.ll) THEN 
               zeile = line (indxb:ll) 
               i = ll-indxb + 1 
            ELSE 
               zeile = ' ' 
               i = 0 
            ENDIF 
            CALL struc_mole_header (zeile, i, .false., lcontent) 
            IF (ier_num.ne.0) return 
         ELSEIF (str_comp (befehl, 'atoms', 2, lbef, 5) ) THEN 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF(ianz==0)   THEN   ! Pre 5.17.2 style, no params
               at_ianz = 4       ! At least x,y,z,Biso
               at_param(1) = 'X'
               at_param(2) = 'Y'
               at_param(3) = 'Z'
               at_param(4) = 'BISO'
               at_param(5:) = ' '
            ELSE
               DO i=1,ianz
                  CALL do_cap(cpara(i))
                  at_param(i) = cpara(i)(1:MIN(LEN(at_param),lpara(i)))
               ENDDO
               at_ianz = ianz
            ENDIF
         ELSE 
            ier_num = - 89 
            ier_typ = ER_APPL 
            RETURN
         ENDIF 
         IF (ier_num.ne.0) THEN 
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
         IF(ier_num/=0) RETURN
      ELSE 
!                                                                       
!     Read old header                                                   
!                                                                       
         ier_num = - 47 
         ier_typ = ER_APPL 
         ier_msg(1) = 'Error in a non-keyword structure file'
         ier_msg(2) = 'Check if the 1st line starts with ''title'''
         READ (ist, 2010, end = 999, err = 999) line 
         lp = len_str(line) 
         call tab2blank(line, lp)
         CALL get_params (line, ianz, cpara, lpara, maxw, lp) 
         cr_spcgr = cpara (1) (1:lpara(1))
         cr_set   = 'abc'
         ianz = ianz - 1 
         IF (ianz.eq.1) THEN 
            cpara (1) = cpara (2) 
            lpara (1) = lpara (2) 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         ENDIF 
         IF (ier_num.eq.0) THEN 
            ier_num = - 48 
            ier_typ = ER_APPL 
            READ (ist, *, end = 999, err = 999) (cr_a0 (i), i = 1, 3),  &
            (cr_win (i), i = 1, 3)                                      
            CALL no_error 
            CALL spcgr_no (ianz, maxw, werte) 
            IF(ier_num/=0) RETURN
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
SUBROUTINE struc_read_atoms (NMAX, MAXSCAT, cr_natoms, cr_nscat,        &
      cr_dw, cr_occ, cr_at_lis, cr_pos, cr_iscat, cr_mole, cr_surf,     &
      cr_magn, cr_prop, cr_dim, cr_magnetic,                            &
      as_natoms, as_at_lis, as_dw, as_occ, as_pos, as_iscat, as_prop,   &
      AT_MAXP, at_ianz, at_param)                      
!-                                                                      
!           This subroutine reads the list of atoms into the            
!       crystal array                                                   
!+                                                                      
USE discus_allocate_appl_mod , ONLY: alloc_molecule
USE molecule_mod 
USE prop_para_mod
USE spcgr_apply
use atom_line_mod
!
use blanks_mod
USE lib_errlist_func
USE lib_length
USE precision_mod
USE str_comp_mod
USE string_convert_mod
!
IMPLICIT none 
!                                                                       
INTEGER                                   , INTENT(IN)    :: NMAX 
INTEGER                                   , INTENT(IN)    :: MAXSCAT 
INTEGER                                   , INTENT(INOUT) :: cr_natoms
INTEGER                                   , INTENT(INOUT) :: cr_nscat 
REAL(kind=PREC_DP), DIMENSION(0:MAXSCAT)  , INTENT(INOUT) :: cr_dw       ! (0:MAXSCAT) 
REAL(kind=PREC_DP), DIMENSION(0:MAXSCAT)  , INTENT(INOUT) :: cr_occ      ! (0:MAXSCAT) 
CHARACTER(LEN=4)  , DIMENSION(0:MAXSCAT)  , INTENT(INOUT) :: cr_at_lis   ! (0:MAXSCAT) 
REAL(kind=PREC_DP), DIMENSION(1:3,1:NMAX) , INTENT(INOUT) :: cr_pos
INTEGER           , DIMENSION(1:NMAX),      INTENT(INOUT) :: cr_iscat
INTEGER           , DIMENSION(1:NMAX),      INTENT(INOUT) :: cr_mole 
INTEGER           , DIMENSION(0:3,1:NMAX) , INTENT(INOUT) :: cr_surf
REAL(kind=PREC_DP), DIMENSION(0:3,1:NMAX) , INTENT(INOUT) :: cr_magn
INTEGER           , DIMENSION(1:NMAX),      INTENT(INOUT) :: cr_prop
REAL(kind=PREC_DP), DIMENSION(3, 2)       , INTENT(INOUT) :: cr_dim      ! (3, 2) 
LOGICAL                                   , INTENT(INOUT) :: cr_magnetic
INTEGER                                   , INTENT(INOUT) :: as_natoms 
CHARACTER(LEN=4)  , DIMENSION(0:MAXSCAT)  , INTENT(INOUT) :: as_at_lis   ! (0:MAXSCAT) 
REAL(kind=PREC_DP), DIMENSION(0:MAXSCAT)  , INTENT(INOUT) :: as_dw       ! (0:MAXSCAT) 
REAL(kind=PREC_DP), DIMENSION(0:MAXSCAT)  , INTENT(INOUT) :: as_occ      ! (0:MAXSCAT) 
REAL(kind=PREC_DP), DIMENSION(3,1:MAXSCAT), INTENT(INOUT) :: as_pos      ! (3, MAXSCAT) 
INTEGER           , DIMENSION(1:MAXSCAT)  , INTENT(INOUT) :: as_iscat    ! (MAXSCAT) 
INTEGER           , DIMENSION(1:MAXSCAT)  , INTENT(INOUT) :: as_prop     ! (MAXSCAT) 
INTEGER                                   , INTENT(IN)    :: AT_MAXP
INTEGER                                   , INTENT(IN )   :: at_ianz
CHARACTER(LEN=8)  , DIMENSION(AT_MAXP)    , INTENT(IN )   :: at_param
!                                                                       
INTEGER , PARAMETER :: ist  = 7
INTEGER , PARAMETER :: maxw = 16! SHOULD READ : MAX(7, AT_MAXP)
!                                                                       
CHARACTER(LEN=10)   :: befehl 
CHARACTER(LEN=PREC_STRING) ::  line, zeile 
INTEGER             :: i, j, ibl, lbef 
INTEGER             :: iatom
INTEGER             :: iimole
INTEGER             :: lline 
INTEGER             :: n_gene
INTEGER             :: n_symm
INTEGER             :: n_mole
INTEGER             :: n_type
INTEGER             :: n_atom
INTEGER             :: n_mole_old = 0    ! previous number of molecules in structure
LOGICAL             :: need_alloc = .false.
LOGICAL             :: lcontent
LOGICAL, SAVE       :: at_init = .TRUE.
REAL, PARAMETER     :: eps = 1e-6
REAL(KIND=PREC_DP), DIMENSION(maxw) :: werte !(maxw)
REAL                :: dw1 , occ1 = 1
!                                                                       
!                                                                       
      lcontent = .false.
      at_init = .TRUE.
      IF(cr_natoms == 0) THEN 
         n_mole_old = 0  ! For empty crystal reset number of old molecules
      ELSE
         n_mole_old = mole_num_mole
      ENDIF
werte = 0.0
 1000 CONTINUE 
      ier_num = - 49 
      ier_typ = ER_APPL 
      line = ' ' 
      READ (ist, 2000, end = 2, err = 999) line 
      lline = len_str (line) 
      call tab2blank(line, lline)
      IF (line.ne.' '.and.line (1:1) .ne.'#'.and. line(1:1) /= '!' .AND. &
          line.ne.char (13) )  THEN
         ibl = index (line (1:lline) , ' ') + 1 
         lbef = 10 
         befehl = ' ' 
         ibl = index (line (1:lline) , ' ') 
         IF (ibl.eq.0) THEN 
            ibl = lline+1 
         ENDIF 
         lbef = min (ibl - 1, lbef) 
         befehl = line (1:lbef) 
         IF (str_comp (befehl, 'molecule', 4, lbef, 8) .OR. &
             str_comp (befehl, 'domain',   4, lbef, 6) .OR. &
             str_comp (befehl, 'object',   4, lbef, 6)     ) THEN
!                                                                       
!     ------Start/End of a molecule                                     
!                                                                       
            CALL no_error 
            IF (ibl.le.lline) THEN 
               zeile = line (ibl:lline) 
               i = lline-ibl + 1 
            ELSE 
               zeile = ' ' 
               i = 0 
            ENDIF 
            CALL struc_mole_header (zeile, i, .false., lcontent) 
            IF (ier_num.ne.0) return 
         ELSE 
!                                                                       
!        --Make sure we have enough space for molecule atoms
!                                                                       
            need_alloc = .false.
            n_gene = MAX( 1, MOLE_MAX_GENE)
            n_symm = MAX( 1, MOLE_MAX_SYMM)
            n_mole =         MOLE_MAX_MOLE
            n_type =         MOLE_MAX_TYPE
            n_atom =         MOLE_MAX_ATOM
            IF ((mole_off(mole_num_mole)+mole_len(mole_num_mole)) >= MOLE_MAX_ATOM ) THEN
               n_atom = (mole_off(mole_num_mole)+mole_len(mole_num_mole)) + 200
               need_alloc = .true.
            ENDIF
            IF ( need_alloc ) THEN
               call alloc_molecule(n_gene, n_symm, n_mole, n_type, n_atom)
            ENDIF
            cr_natoms = cr_natoms + 1 
            cr_surf(:,cr_natoms) = 0
            cr_magn(:,cr_natoms) = 0.0
            call get_atom_werte(cr_natoms, MAXW, werte)
            IF (ier_num.ne.0.and.ier_num.ne. - 49) THEN 
               GOTO 999 
            ENDIF 
            IF (cr_natoms.eq.nmax+1) THEN 
!                                                                       
!     --------Too many atoms in the structure file                      
!                                                                       
               ier_num = - 10 
               ier_typ = ER_APPL 
               RETURN
            ENDIF 
            i = cr_natoms 
            DO j = 1, 3 
            cr_pos (j, i) = werte (j) 
            cr_dim (j, 1) = min (cr_dim (j, 1), cr_pos (j, i) ) 
            cr_dim (j, 2) = max (cr_dim (j, 2), cr_pos (j, i) ) 
            ENDDO 
            dw1 = werte (4) 
            occ1 = werte(8)                             ! WORK OCC
            cr_surf(0:3,i) = NINT(werte(9:12))          ! copy surface
            cr_magn(0:3,i) = werte(13:16)               ! MAGNETIC_WORK
            IF(MAXVAL(ABS(werte(13:16)))>0.0) cr_magnetic = .TRUE.
            cr_prop (i) = nint (werte (5) ) 
      IF (line (1:ibl) .ne.'    ') THEN 
               ibl = ibl - 1 
               CALL do_cap (line (1:ibl) ) 
               DO j = 0, cr_nscat 
                  IF (line (1:ibl) .eq.cr_at_lis (j) .and. &
                      ABS(dw1-cr_dw(j)).lt.eps       .AND. &
                      ABS(occ1-cr_occ(j))<eps             ) THEN
                     cr_iscat (i) = j 
                     GOTO 11 
                  ENDIF 
               ENDDO 
               IF (cr_nscat.eq.MAXSCAT) THEN 
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
               cr_occ(cr_nscat) = occ1    ! WORK OCC
!                                                                       
               IF (0.0.le.cr_pos (1, i) .and.cr_pos (1, i)              &
               .lt.1.and.0.0.le.cr_pos (2, i) .and.cr_pos (2, i)        &
               .lt.1.and.0.0.le.cr_pos (3, i) .and.cr_pos (3, i) .lt.1) &
               THEN                                                     
                  as_natoms = as_natoms + 1 
!                 as_at_lis (cr_nscat) = cr_at_lis (cr_nscat) 
!                 as_iscat (as_natoms) = cr_iscat (i) 
!                 as_dw (as_natoms) = cr_dw (cr_nscat) 
!                 as_occ(as_natoms) = cr_occ(cr_nscat) 
!                 DO j = 1, 3 
!                 as_pos (j, as_natoms) = cr_pos (j, i) 
!                 ENDDO 
!                 as_prop (as_natoms) = cr_prop (i) 
               ENDIF 
   11          CONTINUE 
!                                                                       
!     --------If we are reading a molecule insert atom into current     
!                                                                       
               IF (mole_l_on) THEN 
                  CALL mole_insert_current (cr_natoms, mole_num_curr) 
                  IF (ier_num.lt.0.and.ier_num.ne. - 49) THEN 
                     GOTO 999 
                  ENDIF 
                  cr_prop(cr_natoms) = IBSET(cr_prop(cr_natoms),PROP_MOLECULE)
                  cr_mole(cr_natoms) = mole_num_curr
               ELSE             ! No molecule header, but explicit info on line
                  IF(NINT(werte(6))>0 .AND. NINT(werte(7))>0) THEN
                     iimole = NINT(werte(6)) + n_mole_old
                     CALL mole_insert_explicit(cr_natoms, iimole        , NINT(werte(7))) 
                     cr_prop(cr_natoms) = IBSET(cr_prop(cr_natoms),PROP_MOLECULE)
                     cr_mole(cr_natoms) = NINT(werte(6)) + n_mole_old
                  ENDIF 
               ENDIF 
            ENDIF 
         ENDIF 
      ENDIF 
      GOTO 1000 
!                                                                       
    2 CONTINUE 
!                                                                       
  999 CONTINUE 
!
!     If a molecule containted a "molecule atoms" instruction, we need to
!     set the molecule flag
!
      IF(lcontent) THEN
         DO i = 1, mole_num_mole
            DO j = 1, mole_len (i)
               iatom          = mole_cont (mole_off (i) + j)
               cr_prop(iatom) = ibset(cr_prop(iatom),PROP_MOLECULE)
               cr_mole(iatom) = i
            ENDDO
         ENDDO
      ENDIF 
      IF (ier_num.eq. - 49) THEN 
         CALL no_error 
      ENDIF 
!
!     cr_surf(:,:) = 0    ! Currently no surface save nor read
!
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
USE discus_config_mod 
USE crystal_mod 
USE spcgr_mod 
USE spcgr_apply, ONLY : spcgr_get_setting
!
USE ber_params_mod
USE errlist_mod
USE lib_errlist_func
USE lib_length
USE precision_mod
!
IMPLICIT none 
!                                                                       
INTEGER              , INTENT(INOUT) :: ianz
INTEGER              , INTENT(IN)    :: MAXW
REAL(KIND=PREC_DP), DIMENSION(MAXW), INTENT(IN)    :: werte
!                                                                       
CHARACTER(LEN=   3), DIMENSION(6) :: setting 
CHARACTER(LEN=PREC_STRING), DIMENSION(1) :: cpara 
INTEGER            , DIMENSION(1) :: lpara 
REAL(KIND=PREC_DP) , DIMENSION(1) :: rpara 
!
INTEGER :: ii, i 
INTEGER :: j
!                                                                       
!
DATA setting /'abc', 'bac', 'cab', 'cba', 'bca', 'acb'/
!                                                                       
CALL no_error 
ii = 1 
IF (ianz.eq.1) THEN 
   ii = nint (werte (1) ) 
ENDIF 
!                                                                       
!     Distinguish between trigonal and rhombohedral settings            
!                                                                       
IF(cr_spcgr(1:1) .eq.'R') THEN 
   IF (cr_a0(1) == cr_a0(2) .AND. cr_win(1) == 90. .AND. &
       cr_win(2) == 90.0 .AND. cr_win(3) == 120.0       ) THEN                    
      ii = 1 
   ELSEIF(cr_a0(1) == cr_a0(2) .AND. cr_a0(1) == cr_a0(3) .AND.  &
          cr_win(1) == cr_win(2) .AND. cr_win(1) == cr_win(3)  ) THEN                                                           
      ii = 2 
   ELSE 
      ier_num = - 7 
      ier_typ = ER_APPL 
   ENDIF 
ENDIF 
!                                                                       
IF (ii.lt.1.or.2.lt.ii) THEN 
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
!  Test regular tabulated space groups
!
DO i = 1, SPCGR_MAX 
   IF(cr_spcgr == spcgr_name(i) ) THEN 
      cr_spcgrno = spcgr_num (i, ii) 
      cr_syst = spcgr_syst (cr_spcgrno) 
      CALL no_error 
      GOTO 10 
   ENDIF 
ENDDO 
!
! Test for alternative settings in orthorhombic space groups
!                                                                       
DO i = 16,74
   DO j = 1,6
      IF (cr_spcgr.eq.spcgr_name_set (i,j) ) THEN 
         cr_spcgrno   = i                       ! Here we can use number i as spcgr number instead of spcgr_num (i, ii) 
         cr_syst      = spcgr_syst (cr_spcgrno) ! Still works
         cr_spcgr     = spcgr_name_set(i,1)     ! Use standard name
         cr_spcgr_set = spcgr_name_set(i,j)     ! Alternative setting
         cr_set       = setting(j)
         cr_iset      = j
         CALL no_error 
         GOTO 10 
      ENDIF 
   ENDDO 
ENDDO 
!                                                                       
CALL no_error 
cpara(1) = cr_spcgr 
lpara = len_str (cpara(1)) 
CALL ber_params (1, cpara, lpara, rpara, 1) 
IF (ier_num.eq.0) THEN 
   cr_spcgrno = spcgr_num(nint (rpara(1)) , ii)
   cr_syst = spcgr_syst (cr_spcgrno) 
   cr_spcgr = spcgr_name (cr_spcgrno) 
!
!     Ensure that space groups given as number without origin choice 2 
!     are set properly
!
   IF((276<=cr_spcgrno .AND. cr_spcgrno<=293) .OR. &
      (301<=cr_spcgrno .AND. cr_spcgrno<=306)     ) THEN
      spcgr_para = 2
      spcgr_ianz = 1
   ENDIF
ELSE 
   ier_num = - 7 
   ier_typ = ER_APPL 
ENDIF 
!                                                                       
10 CONTINUE 
!                                                                       
IF (ier_num == 0) THEN 
   ier_num = - 14 
   ier_typ = ER_APPL 
   IF (cr_syst == 1) THEN 
      CALL no_error 
   ELSEIF(cr_syst == 2) THEN 
      IF(cr_win(1) == 90.0 .AND. cr_win(3) == 90.0) THEN 
         CALL no_error 
      ENDIF 
   ELSEIF(cr_syst == 3) THEN 
      IF(cr_win(1) == 90.0 .AND. cr_win(2) == 90.0) THEN 
         CALL no_error 
      ENDIF 
   ELSEIF(cr_syst == 4) THEN 
      IF(cr_win(1) == 90.0 .AND. cr_win(2) == 90.0 .AND. cr_win(3) == 90.0) THEN                                           
         CALL no_error 
      ENDIF 
   ELSEIF(cr_syst == 5) THEN 
      IF(cr_a0(1) == cr_a0(2) .AND. cr_win(1) == 90.0 .AND. &
         cr_win(2) == 90.0    .AND. cr_win(3) == 90.0      )   THEN                                                        
         CALL no_error 
      ENDIF 
   ELSEIF(cr_syst == 6 .OR. cr_syst == 8) THEN 
      IF(cr_a0(1) == cr_a0(2) .AND. cr_win(1) == 90.0 .AND. &
         cr_win(2) == 90.0    .AND. cr_win(3) == 120.0     )  THEN                                                        
         CALL no_error 
      ENDIF 
   ELSEIF(cr_syst == 7) THEN 
      IF(cr_a0(1) == cr_a0(2)   .AND. cr_a0(1) == cr_a0(3) .AND.   &
         cr_win(1) == cr_win(2) .AND. cr_win(2) == cr_win(3)    ) THEN                                                   
         CALL no_error 
      ENDIF 
   ELSEIF(cr_syst == 9) THEN 
      IF (cr_a0(1) == cr_a0(2) .AND. cr_a0(1) == cr_a0(3)  .AND.  &
          cr_win(1) == 90.0    .AND. cr_win(2) == 90.0     .AND.  &
          cr_win(3) == 90.0                                     ) THEN                                          
         CALL no_error 
      ENDIF 
   ENDIF 
ENDIF
!
CALL spcgr_get_setting    ! Determine space group setting
IF(ier_num/=0) RETURN
IF(cr_syst/=4 .AND. cr_iset/=1) THEN !non-orthorhombic and nonstandard setting
   ier_num = -160
   ier_typ = ER_APPL
   ier_msg(1) = 'Non-standard settings are implemented only'
   ier_msg(2) = 'for orthorhombic space groups '
   RETURN
ENDIF
!
END SUBROUTINE spcgr_no                       
!********************************************************************** 
!********************************************************************** 
SUBROUTINE rese_cr 
!                                                                       
!     resets the crystal structure to empty                             
!                                                                       
USE discus_config_mod 
USE discus_allocate_appl_mod
USE conn_mod
USE crystal_mod 
USE gen_add_mod 
USE molecule_mod 
USE sym_add_mod 
USE discus_save_mod 
USE discus_plot_mod
USE precision_mod
!
IMPLICIT none 
!                                                                       
INTEGER, PARAMETER  :: code_res   = -2
!
CHARACTER(LEN=PREC_STRING) :: zeile
INTEGER             :: lp
!                                                                       
      INTEGER i 
!
CALL alloc_crystal(1,1)                                                                       
!
      cr_natoms = 0 
      as_natoms = 0 
      cr_ncatoms = 1 
      cr_ncreal  = 1 
      cr_nscat = 0 
      cr_icc       = 1
      cr_cartesian = .false. 
cr_magnetic  = .FALSE.
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
      cr_dw(:)     = 0.0
      cr_occ(:)    = 1.0
      as_dw(:)     = 0.0
      as_occ(:)    = 1.0
!
      DO i = 1, 3 
      cr_dim (i, 1) = 0.0D0 
      cr_dim (i, 2) = 0.0D0
      ENDDO 
!
cr_amount = 0
cr_u2aver = 0.0E0
cr_mass   = 0   ! Crystal mass in u
cr_nreal  = 0.0
      as_pos(:,:)  = 0
      cr_pos(:,:)  = 0
      cr_prop(:)   = 0
      cr_iscat(:)  = 0
      cr_niscat(:)  = 0
      cr_mole(:)   = 0
      cr_surf(:,:) = 0
      cr_magn(:,:) = 0.0
!                                                                       
!     DO i = 0, MOLE_MAX_MOLE 
      mole_len (:) = 0 
      mole_off (:) = 0 
      mole_type(:) = 0 
      mole_char(:) = 0 
      mole_file(:) = ' '
      mole_cont(:) = 0 
      mole_dens(:) = 0.0
      mole_biso(:) = 0.0
      mole_clin(:) = 0.0
      mole_cqua(:) = 0.0
      mole_fuzzy(:) = 0.0
      mole_gene (:,:,:) = 0.0
      mole_symm (:,:,:) = 0.0
!     ENDDO 
      mole_l_on = .false. 
      mole_num_mole = 0 
      mole_num_curr = 0 
      mole_num_act = 0 
      mole_num_type = 0 
      mole_num_unit = 0 
      mole_num_atom = 0 
      mole_num_acur = 0 
      mole_gene_n = 0 
      mole_symm_n = 0 
!                                                                       
      sym_add_n = 0 
      gen_add_n = 0 
!                                                                       
      sav_r_ncell = .false. 
!
      pl_poly_n = 0
      pl_poly_c(:) = .FALSE.
      pl_poly_o(:) = .FALSE.
      pl_poly_dmin = 0.0
      pl_poly_dmax = 0.0
      pl_poly_nmin = 0
      pl_poly_dmax = 0
      pl_poly_face = .TRUE.
      pl_poly_hue  = .FALSE.
      pl_poly_col  = 'auto'
!
      pl_bond(:,:)       = .FALSE.
      pl_bond_len(:,:,:) = 0.0
      pl_bond_rad(  :,:) = 0.0
      pl_bond_col(:,:,:) = 0.0
!
zeile = ' '
lp    = 1
CALL conn_do_set(code_res,zeile, lp)    ! Connectivity
!                                                                       
      END SUBROUTINE rese_cr                        
!
!*****7**************************************************************** 
!
SUBROUTINE import_test(mode, strufile, outfile)
!
! Tests if the ending of a file corresponds to a known format, tries to 
! import this file
!
USE lib_errlist_func
USE errlist_mod
INTEGER         , INTENT(IN)  :: mode
CHARACTER(LEN=*), INTENT(IN)  :: strufile
CHARACTER(LEN=*), INTENT(OUT) ::  outfile
!
!INTEGER, PARAMETER :: MD_CELL = 0
INTEGER, PARAMETER :: MD_STRU = 1
!
CHARACTER(LEN=LEN_TRIM(strufile)+8) :: line
INTEGER :: length
INTEGER :: laenge
!
length = LEN_TRIM(strufile)
!
IF(strufile(length-3:length) == '.cif' .OR. strufile(length-3:length) == '.CIF') THEN
   line = 'cif, '//strufile
   laenge = 5 + length
   CALL do_import(line, laenge)
   IF(ier_num == 0) THEN
      outfile = strufile(1:length-3) // 'stru'
   ENDIF
ELSEIF(strufile(length-3:length) == '.txt' .OR. strufile(length-3:length) == '.TXT') THEN
   line = 'cmaker, '//strufile
   laenge = 8 + length
   CALL do_import(line, laenge)
   IF(ier_num == 0) THEN
      outfile = strufile(1:length-3) // 'stru'
   ENDIF
ELSEIF(strufile(length-4:length) == '.cssr' .OR. strufile(length-4:length) == '.CSSR') THEN
   IF(mode==MD_STRU ) THEN
      line = 'rmc, '//strufile
      laenge = 5 + length
      CALL do_import(line, laenge)
      IF(ier_num == 0) THEN
         outfile = strufile(1:length-4) // 'stru'
      ENDIF
   ELSE
      ier_num = -140
      ier_typ = ER_APPL
      ier_msg(1) = 'RMCprofile files usually contain more than'
      ier_msg(2) = 'one unit cell'
      ier_msg(3) = 'import first and convert the unit cell size'
   ENDIF
ELSE
    CALL no_error
    outfile = strufile
ENDIF
!
END SUBROUTINE import_test
!*****7**************************************************************** 
!
SUBROUTINE do_import(zeile, lp) 
!-                                                                      
!     imports a file into discus.cell format                            
!+                                                                      
USE crystal_mod
USE prop_para_func
USE discus_save_mod
USE read_internal_mod
USE save_menu, ONLY: save_internal, save_store_setting, &
                     save_restore_setting, save_default_setting, &
                     save_struc, save_show, save_keyword
USE spcgr_apply
!USE trafo_mod
!
USE build_name_mod
USE get_params_mod
USE lib_errlist_func
USE precision_mod
USE str_comp_mod
USE take_param_mod
!                                                                       
IMPLICIT none 
!                                                                       
CHARACTER(LEN=*), INTENT(INOUT) :: zeile 
INTEGER         , INTENT(INOUT) :: lp 
!                                                                       
INTEGER, PARAMETER :: MAXW = 5 
!                                                                       
CHARACTER(LEN=MAX(PREC_STRING,LEN(ZEILE))) ::  ofile
CHARACTER(LEN=MAX(PREC_STRING,LEN(ZEILE))) ::  line
CHARACTER(LEN=MAX(PREC_STRING,LEN(zeile))), DIMENSION(MAXW) ::  cpara !(MAXW) 
INTEGER            , DIMENSION(MAXW) ::  lpara !(MAXW) 
REAL(KIND=PREC_DP) , DIMENSION(MAXW) ::  werte !(MAXW) 
INTEGER :: ianz 
INTEGER :: length
!
CHARACTER(LEN= 200) :: hostfile  ! original structure file name
!
INTEGER, PARAMETER :: NOPTIONAL = 4
INTEGER, PARAMETER :: O_METRIC  = 1
INTEGER, PARAMETER :: O_SPACE   = 2
INTEGER, PARAMETER :: O_SORT    = 3
INTEGER, PARAMETER :: O_ATOM    = 4     ! For LAMMPS
CHARACTER(LEN=   6)       , DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=MAX(PREC_STRING,LEN(zeile))), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para is present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 0 ! Number of values to calculate 
!
DATA oname  / 'metric', 'space', 'sort', 'atom'   /
DATA loname /  6,        5     ,  4    ,  4       /
!
REAL, DIMENSION(1:3) :: host_a0
REAL, DIMENSION(1:3) :: host_win
REAL, DIMENSION(4,4) :: host_tran_fi
INTEGER              :: host_syst
INTEGER              :: host_spcgrno
CHARACTER(LEN=16)    :: host_spcgr
CHARACTER(LEN= 3)    :: host_set   = 'abc'
INTEGER              :: host_iset  =  1
CHARACTER(LEN=16)    :: host_spcgr_set = 'P1'
REAL, DIMENSION(4)   :: posit4 ! atom position
REAL, DIMENSION(4)   :: uvw4   ! atom position
INTEGER              :: j
LOGICAL              :: lperiod
!
LOGICAL :: lout = .FALSE.
!
!
opara  =  (/ 'guest ', 'P1    ', 'discus', 'atom  ' /)   ! Always provide fresh default values
lopara =  (/  6,        6      ,  6      ,  6       /)
owerte =  (/  0.0,      0.0    ,  0.0    ,  0.0     /)
!                                                                       
CALL get_params (zeile, ianz, cpara, lpara, MAXW, lp) 
IF (ier_num.ne.0) THEN 
   RETURN
ENDIF 
!                                                                       
CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
IF (ier_num.ne.0) THEN 
   RETURN
ENDIF 
!
lperiod=.FALSE.
!
IF (ianz.ge.1) THEN 
   IF (str_comp (cpara (1) , 'shelx', 2, lpara (1) , 5) ) THEN 
      IF (ianz >= 2) THEN 
         CALL del_params (1, ianz, cpara, lpara, maxw) 
         IF (ier_num.ne.0) return 
         CALL ins2discus (ianz, cpara, lpara, MAXW, ofile) 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
   ELSEIF (str_comp (cpara (1) , 'cif', 2, lpara (1) , 3) ) THEN 
      IF (ianz >= 2) THEN 
         CALL del_params (1, ianz, cpara, lpara, maxw) 
         IF (ier_num.ne.0) return 
         CALL cif2discus (ianz, cpara, lpara, MAXW, ofile) 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
   ELSEIF (str_comp (cpara (1) , 'cmaker', 2, lpara (1) , 6) ) THEN 
      IF (ianz >= 2) THEN 
         CALL del_params (1, ianz, cpara, lpara, maxw) 
         IF (ier_num.ne.0) return 
         CALL cmaker2discus (ianz, cpara, lpara, MAXW, ofile) 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
   elseif(str_comp(cpara(1), 'lammps', 2, lpara(1) , 6) ) THEN 
      if (ianz >= 2) then 
         call del_params(1, ianz, cpara, lpara, maxw) 
         if (ier_num /= 0) return 
         call lammps2discus(ianz, cpara, lpara, MAXW,           &
                            NOPTIONAL, opara, lopara, lpresent, O_ATOM, ofile) 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
   ELSEIF (str_comp (cpara (1) , 'rmcprofile', 2, lpara (1) , 10) ) THEN 
      IF(lpresent(O_SORT)) THEN
         IF(str_comp(opara(O_SORT), 'discus', 6, lopara(O_SORT), 6)) THEN
            lperiod=.TRUE.
         ELSE
            lperiod=.FALSE.
         ENDIF
      ELSE
         lperiod = .FALSE.
      ENDIF
      IF (ianz >= 2) THEN 
         CALL del_params (1, ianz, cpara, lpara, maxw) 
         IF (ier_num.ne.0) return 
         CALL rmcprofile2discus (ianz, cpara, lpara, MAXW, ofile, lperiod) 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
   ELSEIF(str_comp(cpara(1), 'cell', 2, lpara(1), 4) .OR.          &      
          str_comp(cpara(1), 'stru', 2, lpara(1), 4)     ) THEN
      IF(opara(O_METRIC)=='host' ) THEN 
         IF (ianz >= 2) THEN 
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            IF (ier_num.ne.0) RETURN 
            CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
            IF (ier_num.ne.0) THEN 
               RETURN
            ENDIF 
            ofile = cpara(1)(1:lpara(1))
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ELSE 
         ier_num = -167
         ier_typ = ER_COMM 
         ier_msg(1) = 'Provide the optional parameter'
      ENDIF 
   ELSE 
      ier_num = - 86 
      ier_typ = ER_APPL 
   ENDIF 
ELSE 
   ier_num = - 6 
   ier_typ = ER_COMM 
ENDIF 
IF(ier_num==0) THEN
  IF(opara(O_METRIC)=='host')  THEN  ! Transform into lattice parameter of host
!
   CALL save_store_setting             ! Backup user "save" setting
   CALL save_default_setting           ! Default to full saving
   line       = 'ignore, all'          ! Ignore all properties
   length     = 11
   CALL property_select(line, length, sav_sel_prop)
   line       = 'ignore, all'          ! Ignore all properties for global as well
   length     = 11
   CALL property_select(line, length,  cr_sel_prop)
   hostfile = 'internal.import'
   CALL save_internal(hostfile)        !     thus this file name is unique
!
   host_a0        = cr_a0              ! Store host lattice parameters
   host_win       = cr_win             ! same for angles
   host_tran_fi   = cr_tran_fi
   host_syst      = cr_syst
   host_spcgr     = cr_spcgr
   host_spcgrno   = cr_spcgrno
   host_set       = cr_set
   host_iset      = cr_iset
   host_spcgr_set = cr_spcgr_set
!
   CALL rese_cr 
   CALL do_readstru_disk(ofile, .FALSE.)
   CALL setup_lattice (cr_a0, cr_ar, cr_eps, cr_gten, cr_reps, &
        cr_rten, cr_win, cr_wrez, cr_v, cr_vr, lout, cr_gmat,  &
        cr_fmat, cr_cartesian,                                 &
        cr_tran_g, cr_tran_gi, cr_tran_f, cr_tran_fi)
   CALL get_symmetry_matrices 
!
! Transform molecule geometry into host geometry
   DO j=1, cr_natoms
      posit4(1:3) = cr_pos(1:3,j)
      posit4(4)   = 1.0
!     CALL trans(posit4,cr_tran_f ,uvw4, 4)       ! Transform mol to cartesian
!     CALL trans(uvw4  ,host_tran_fi,posit4, 4)   ! Transform cartesian to host_crystal
      uvw4   = matmul(cr_tran_f, posit4)          ! Transform mol to cartesian
      posit4 = matmul(host_tran_fi, uvw4)         ! Transform cartesian to host_crysta
      cr_pos(1:3,j) = posit4(1:3)
   ENDDO
!
   cr_a0(:)  = host_a0
   cr_win(:) = host_win
!
   IF(opara(O_SPACE)=='host') THEN
      cr_syst      = host_syst
      cr_spcgr     = host_spcgr
      cr_spcgrno   = host_spcgrno
      cr_spcgr_set = host_spcgr_set
      cr_iset      = host_iset
   ELSE
      cr_syst      = CR_TRICLINIC
      cr_spcgr     = 'P1'
      cr_spcgrno   = 1
      cr_spcgr_set = 'abc'
      cr_iset      = 1
   ENDIF
!
   CALL save_default_setting           ! Default to full saving
   sav_latom(0:cr_nscat) = .TRUE.
   sav_w_scat  = .TRUE.
   sav_w_adp   = .FALSE.
   sav_w_occ   = .FALSE.
   sav_w_surf  = .TRUE.
   sav_w_magn  = .FALSE.  ! MAGNETIC_WORK
   sav_r_ncell = .FALSE.
   sav_w_ncell = .FALSE.
   sav_w_gene  = .FALSE.
   sav_w_symm  = .FALSE.
   sav_w_mole  = .TRUE.
   sav_w_obje  = .TRUE.
   sav_w_doma  = .TRUE.
   sav_w_prop  = .TRUE.
   CALL save_keyword(ofile)
!
   CALL save_restore_setting
   CALL no_error
   CALL readstru_internal( hostfile)   ! Read  core file
   CALL errlist_restore                ! Restore error status

  ENDIF
ENDIF
!                                                                       
END SUBROUTINE do_import                      
!*****7**************************************************************** 
SUBROUTINE ins2discus (ianz, cpara, lpara, MAXW, ofile) 
!-                                                                      
!     converts a SHELXL "ins" or "res" file to DISCUS                   
!+                                                                      
USE ber_params_mod
USE blanks_mod
USE build_name_mod
USE get_params_mod
USE lib_length
USE wink_mod
USE precision_mod
USE support_mod
!
IMPLICIT none 
!                                                                       
!                                                                       
INTEGER                             , INTENT(INOUT) :: ianz 
INTEGER                             , INTENT(IN)    :: MAXW 
CHARACTER (LEN= * ), DIMENSION(MAXW), INTENT(INOUT) :: cpara ! (MAXW) 
INTEGER            , DIMENSION(MAXW), INTENT(INOUT) :: lpara ! (MAXW) 
CHARACTER(LEN=*)                    , INTENT(OUT)   :: ofile 
!                                                                       
INTEGER, PARAMETER :: NFV = 50 
!                                                                       
      REAL(KIND=PREC_DP) :: werte (3) 
!                                                                       
      INTEGER shelx_num 
      PARAMETER (shelx_num = 61) 
      CHARACTER(4) shelx_ign (1:shelx_num) 
      CHARACTER(2) c_atom (20) 
      CHARACTER(4) command 
CHARACTER(LEN=4), DIMENSION(:), ALLOCATABLE :: eadp_names
      CHARACTER(80) line1 
      CHARACTER(80) line2 
      CHARACTER(160) line 
      CHARACTER(LEN=MAX(PREC_STRING,LEN(ofile))) :: infile 
      INTEGER ird, iwr 
      INTEGER i, j, jj 
      INTEGER ix, iy, iz, idot 
      INTEGER ntyp , ntyp_prev
      INTEGER length, length1, length2, lp 
      INTEGER icont 
      INTEGER centering 
      INTEGER ityp 
      INTEGER ifv 
INTEGER   :: MAX_EADP
INTEGER   :: n_eadp
      LOGICAL lmole, lmole_wr 
      LOGICAL lcontinue 
      REAL z, latt (6) 
      REAL xyz (3) 
      REAL uiso, uij (6) 
      REAL gen (3, 4) 
      REAL fv (NFV) 
REAL   , DIMENSION(:), ALLOCATABLE :: eadp_values
!
      INTEGER                               :: iianz      ! Dummy number of parameters
      INTEGER, PARAMETER                    :: MAXP  = 11 ! Dummy number of parameters
      CHARACTER (LEN=MAX(PREC_STRING,LEN(cpara))), DIMENSION(MAXP) :: ccpara     ! Parameter needed for SFAC analysis
      INTEGER             , DIMENSION(MAXP) :: llpara
      REAL(KIND=PREC_DP)  , DIMENSION(MAXP) :: wwerte
!                                                                       
!                                                                       
      DATA shelx_ign / 'ACTA', 'AFIX', 'ANIS', 'BASF', 'BIND', 'BLOC',  &
      'BOND', 'BUMP', 'CGLS', 'CHIV', 'CONF', 'CONN', 'DAMP', 'DANG',   &
      'DEFS', 'DELU', 'DFIX', 'DISP',         'EQIV', 'EXTI', 'EXYZ',   &
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
      lmole    = .false. 
      lmole_wr = .true. 
!
      ntyp      = 0
      ntyp_prev = 0
!                                                                       
      CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
      IF (ier_num.ne.0) THEN 
         RETURN
      ENDIF 
      infile = cpara (1) 
      i = index (infile, '.', .TRUE.) 
      IF (i.eq.0) THEN 
         infile = cpara (1) (1:lpara (1) ) //'.ins' 
         ofile = cpara (1) (1:lpara (1) ) //'.cell' 
      ELSE 
         ofile = cpara (1) (1:i) //'cell' 
      ENDIF 
      ird = 34 
      iwr = 35 
      CALL oeffne (ird, infile, 'old') 
      IF (ier_num.ne.0) THEN 
         RETURN
      ENDIF 
      CALL oeffne (iwr, ofile, 'unknown') 
      IF (ier_num.ne.0) THEN 
         RETURN
      ENDIF 
!
!   Allocate and initialize the EADP array
!
MAX_EADP = 20
ALLOCATE(eadp_names( 1:MAX_EADP))
ALLOCATE(eadp_values(1:MAX_EADP))
!
n_eadp         = 0
eadp_names(:)  = ' '
eadp_values(:) = 0.0
!                                                                       
      lcontinue = .false. 
      READ (ird, 1000, end = 900, err = 900) line1 
      length1 = len_str (line1) 
      IF (length1.gt.0) THEN 
         icont = index (line1, '=') 
         IF (icont.gt.0) THEN 
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
      IF (length.gt.0) THEN 
         command = line (1:4) 
      ELSE 
         command = '    ' 
      ENDIF 
      DO i = 1, shelx_num 
         lcontinue = lcontinue.or.command.eq.shelx_ign (i) 
      ENDDO 
!                                                                       
      DO while (command.ne.'FVAR'.and.command.ne.'MOLE') 
      IF (lcontinue) THEN 
         CONTINUE 
      ELSEIF (command.eq.'TITL') THEN 
         WRITE (iwr, 2000) line (6:length) 
         WRITE (iwr, 2100) 
      ELSEIF (command.eq.'CELL') THEN 
         READ (line (6:length), * ) z, latt 
         WRITE (iwr, 2200) latt 
      ELSEIF (command.eq.'LATT') THEN 
         READ (line (6:length), * ) centering 
         IF (abs (centering) .eq.1) THEN 
            CONTINUE 
         ELSEIF (abs (centering) .eq.2) THEN 
            WRITE (iwr, 2320) 
         ELSEIF (abs (centering) .eq.3) THEN 
            WRITE (iwr, 2330) 
         ELSEIF (abs (centering) .eq.4) THEN 
            WRITE (iwr, 2340) 
            WRITE (iwr, 2341) 
         ELSEIF (abs (centering) .eq.5) THEN 
            WRITE (iwr, 2350) 
         ELSEIF (abs (centering) .eq.6) THEN 
            WRITE (iwr, 2360) 
         ELSEIF (abs (centering) .eq.7) THEN 
            WRITE (iwr, 2370) 
         ENDIF 
         IF (centering.gt.0) THEN 
            WRITE (iwr, 2400) 
         ENDIF 
      ELSEIF (command.eq.'SFAC') THEN 
         j = 5 
         atom_search: DO while (j.lt.length) 
            j = j + 1 
            DO while (j.lt.length.and.line (j:j) .eq.' ') 
               j = j + 1 
            ENDDO 
            IF (j.le.length) THEN 
               ntyp = ntyp + 1 
               c_atom (ntyp) = ' ' 
               i = 0 
               DO while (j.le.length.and.line (j:j) .ne.' ') 
                  i = i + 1 
                  c_atom (ntyp) (i:i) = line (j:j) 
                  j = j + 1 
               ENDDO 
               IF(ntyp == ntyp_prev + 2) THEN
!
!                 This is the second parameter, test if this is a numerical
!                 value. If so only the first parameter is an atom name rest is
!                 the numerical form factor, which we ignore
                  ccpara(1) = c_atom(ntyp)
                  llpara(1) = i
                  iianz     = 1
                  CALL ber_params (iianz, ccpara, llpara, wwerte, MAXP) 
                  IF(ier_num==0) THEN
                     ntyp = ntyp - 1
                     EXIT atom_search
                  ENDIF
                  ier_num = 0
                  ier_typ = ER_NONE
               ENDIF
            ENDIF 
         ENDDO atom_search
!        WRITE (iwr, 2500) (c_atom (i) , ',', i = 1, ntyp - 1) , c_atom &
!        (ntyp)                                                         
      ELSEIF (command.eq.'SYMM') THEN 
         lp = length - 5 
         CALL get_params (line (6:length), ianz, cpara, lpara, maxw, lp) 
         IF (ianz.eq.3) THEN 
            DO i = 1, 3 
            DO jj = 1, 4 
            gen (i, jj) = 0.0 
            ENDDO 
            ix = index (cpara (i) , 'X') 
            IF (ix.gt.0) THEN 
               gen (i, 1) = 1.0 
               IF (ix.gt.1) THEN
                  IF(cpara (i) (ix - 1:ix - 1) .eq.'-') THEN 
                  gen (i, 1) = - 1.0 
                  cpara (i) (ix - 1:ix - 1) = ' ' 
               ELSEIF (cpara (i) (ix - 1:ix - 1) .eq.'+')   &
               THEN                                                     
                  gen (i, 1) = 1.0 
                  cpara (i) (ix - 1:ix - 1) = ' ' 
               ENDIF 
               ENDIF 
               cpara (i) (ix:ix) = ' ' 
            ENDIF 
            iy = index (cpara (i) , 'Y') 
            IF (iy.gt.0) THEN 
               gen (i, 2) = 1.0 
               IF (iy.gt.1) THEN
                  IF(cpara (i) (iy - 1:iy - 1) .eq.'-') THEN 
                  gen (i, 2) = - 1.0 
                  cpara (i) (iy - 1:iy - 1) = ' ' 
               ELSEIF (cpara (i) (iy - 1:iy - 1) .eq.'+')   &
               THEN                                                     
                  gen (i, 2) = 1.0 
                  cpara (i) (iy - 1:iy - 1) = ' ' 
               ENDIF 
               ENDIF 
               cpara (i) (iy:iy) = ' ' 
            ENDIF 
            iz = index (cpara (i) , 'Z') 
            IF (iz.gt.0) THEN 
               gen (i, 3) = 1.0 
               IF (iz.gt.1) THEN
                  IF(cpara (i) (iz - 1:iz - 1) .eq.'-') THEN 
                  gen (i, 3) = - 1.0 
                  cpara (i) (iz - 1:iz - 1) = ' ' 
               ELSEIF (cpara (i) (iz - 1:iz - 1) .eq.'+')   &
               THEN                                                     
                  gen (i, 3) = 1.0 
                  cpara (i) (iz - 1:iz - 1) = ' ' 
               ENDIF 
               ENDIF 
               cpara (i) (iz:iz) = ' ' 
            ENDIF 
            ENDDO 
            DO i = 1, 3 
            idot = index (cpara (i) , '.') 
            IF (idot.eq.0) THEN 
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
      IF (length1.gt.0) THEN 
         icont = index (line1, '=') 
         IF (icont.gt.0) THEN 
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
      IF (length.gt.0) THEN 
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
      atoms: DO while (command.ne.'HKLF') 
      IF (lcontinue) THEN 
         CONTINUE 
      ELSEIF (command.eq.'FVAR') THEN 
         READ (line (6:length), *, end = 800) (fv (i), i = 1, NFV) 
  800    CONTINUE 
      ELSEIF (command.eq.'MOLE') THEN 
         IF (lmole) THEN 
            WRITE (iwr, 4000) 'molecule end' 
            lmole_wr = .true. 
         ELSE 
            lmole = .true. 
            lmole_wr = .true. 
         ENDIF 
         CONTINUE 
      ELSEIF(command=='EADP') THEN
         CONTINUE
      ELSEIF (command.eq.'    ') THEN 
         CONTINUE
      ELSE 
         IF (lmole.and.lmole_wr) THEN 
            WRITE (iwr, 4000) 'molecule' 
            lmole_wr = .false. 
         ENDIF 
!
!        This is an atom, get the parameters from the input line
!
         iianz  = 0
         j      = 5 
         ccpara = ' '
         llpara = 0
         atom_para: DO while (j.lt.length) 
            j = j + 1 
            DO while (j.lt.length.and.line (j:j) .eq.' ') 
               j = j + 1 
            ENDDO 
            IF (j.le.length) THEN 
               iianz = iianz + 1 
               ccpara (iianz) = ' ' 
               i = 0 
               DO while (j.le.length.and.line (j:j) .ne.' ') 
                  i = i + 1 
                  ccpara (iianz) (i:i) = line (j:j) 
                  j = j + 1 
               ENDDO 
               llpara(iianz) = i
            ENDIF 
         ENDDO atom_para
         READ (ccpara(1)(1:llpara(1)),*) ityp
         READ (ccpara(2)(1:llpara(2)),*) xyz(1)
         READ (ccpara(3)(1:llpara(3)),*) xyz(2)
         READ (ccpara(4)(1:llpara(4)),*) xyz(3)
         uij(:) = 0
         DO i=1,iianz - 5
            READ (ccpara(5+i)(1:llpara(5+i)),*) uij(i)
         ENDDO
!        READ (line (6:length), *, end = 850) ityp, xyz, sof, (uij (i), &
!        i = 1, 6)                                                      
! 850    CONTINUE 
         DO i=1,3
            ifv = nint (uij(i)/10.)
            IF(ifv.gt.1) THEN
               uij(i) = (uij(i) - ifv * 10) * fv (ifv)
            ELSEIF (ifv.lt. - 1) THEN
               uij(i) = (abs (uij(i)) + ifv * 10) * (1. - fv(IABS(ifv)))
            ENDIF
         ENDDO
         IF (iianz == 6) THEN 
            uiso = uij (1) 
         ELSE 
            uiso = (uij (1) + uij (2) + uij (3) ) / 3. 
         ENDIF 
         DO i = 1, 3 
         ifv = nint (xyz (i) / 10.) 
         IF (ifv.eq.1) THEN 
            xyz (i) = xyz (i) - 10. 
         ELSEIF (ifv.gt.1) THEN 
            xyz (i) = (xyz (i) - ifv * 10) * fv (ifv) 
         ELSEIF (ifv.lt. - 1) THEN 
            xyz(i) = (ABS(xyz(i)) + ifv * 10) * (1. - fv(IABS(ifv)))
         ENDIF 
         ENDDO 
!         write(iwr,3100) c_atom(ityp),xyz,REAL(ityp)                  
         WRITE (iwr, 3100) c_atom (ityp), xyz, uiso *8.*pi**2
      ENDIF 
!                                                                       
      lcontinue = .false. 
      READ (ird, 1000, end = 900, err = 900) line1 
      length1 = len_str (line1) 
      IF (length1.gt.0) THEN 
         icont = index (line1, '=') 
         IF (icont.gt.0) THEN 
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
         IF (length.gt.0) THEN 
            command = line (1:4) 
         ELSE 
            command = '    ' 
            CYCLE atoms
         ENDIF 
         DO i = 1, shelx_num 
            lcontinue = lcontinue.or.command.eq.shelx_ign (i) 
         ENDDO 
      ENDDO  atoms
!                                                                       
  900 CONTINUE 
!                                                                       
      IF (lmole .AND. .NOT.lmole_wr) THEN   ! Write a final molecule end 
         WRITE (iwr, 4000) 'molecule end' 
      ENDIF 
!                                                                       
      CLOSE (ird) 
      CLOSE (iwr) 
DEALLOCATE(eadp_names)
DEALLOCATE(eadp_values)
!                                                                       
 1000 FORMAT    (a) 
 2000 FORMAT    ('title ',a) 
 2100 FORMAT    ('spcgr P1') 
 2200 FORMAT    ('cell ',5(2x,f9.4,','),2x,f9.4) 
 2320 FORMAT    ('gener  1.0, 0.0, 0.0, 0.5,',                          &
     &                     '    0.0, 1.0, 0.0, 0.5,',                   &
     &                     '    0.0, 0.0, 1.0, 0.5,  1')                
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
 2600 FORMAT    ('gener',3(2X,4(1x,f12.8,',')),' 1.') 
!                                                                       
 3000 FORMAT    ('atoms') 
 3100 FORMAT    (a2,2x,4(2x,f9.5)) 
!                                                                       
 4000 FORMAT    (a) 
!                                                                       
      END SUBROUTINE ins2discus                     
!*****7**************************************************************** 
      SUBROUTINE cmaker2discus (ianz, cpara, lpara, MAXW, ofile) 
!-                                                                      
!     converts a CrystalMaker "xyz" file to DISCUS                   
!+                                                                      
      USE build_name_mod
USE lib_length
USE precision_mod
USE str_comp_mod
USE support_mod
!
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER          , INTENT(INOUT)                    :: ianz 
      INTEGER          , INTENT(IN)                       :: MAXW 
      CHARACTER (LEN=*), DIMENSION(1:MAXW), INTENT(INOUT) :: cpara
      INTEGER          , DIMENSION(1:MAXW), INTENT(INOUT) :: lpara
CHARACTER(LEN=*)                 , INTENT(OUT)   :: ofile 
!                                                                       
!                                                                       
      REAL(KIND=PREC_DP)   , DIMENSION(3) :: werte
!                                                                       
      CHARACTER(LEN=87)     :: line 
      CHARACTER(LEN=MAX(PREC_STRING,LEN(ofile)))   :: infile 
      INTEGER               :: ird, iwr 
      INTEGER               :: i
      INTEGER               :: indx1, indx2
      INTEGER               :: iostatus
      INTEGER               :: natoms
      INTEGER               :: nline
      INTEGER               :: length
      REAL   , DIMENSION(6) :: latt (6) 
!                                                                       
!                                                                       
!     Create input / output file name
!
      CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
      IF (ier_num.ne.0) THEN 
         RETURN
      ENDIF 
      infile = cpara (1) 
      i = index (infile, '.') 
      IF (i.eq.0) THEN 
         infile = cpara (1) (1:lpara (1) ) //'.txt' 
         ofile  = cpara (1) (1:lpara (1) ) //'.cell' 
      ELSE 
         ofile  = cpara (1) (1:i) //'cell' 
      ENDIF 
      ird = 34 
      iwr = 35 
      CALL oeffne (ird, infile, 'old') 
      IF (ier_num.ne.0) THEN 
         RETURN
      ENDIF 
      CALL oeffne (iwr, ofile, 'unknown') 
      IF (ier_num.ne.0) THEN 
         RETURN
      ENDIF 
!                                                                       
      nline     = 1
header: DO
         READ (ird, 1000, iostat=iostatus) line
         IF(iostatus /= 0) THEN
            CLOSE(ird)
            CLOSE(iwr)
            ier_msg(1) = 'Error reading CrystalMaker file'
            WRITE(ier_msg(2),5000) nline
            RETURN
         ENDIF
         length = len_str (line) 
zero_o:  IF (length.gt.0) THEN 
cmd:        IF(str_comp(line(1:4),'Unit', 4, length, 4)) THEN
               READ (ird, 1000, iostat=iostatus) line
               nline = nline + 1
               IF(iostatus /= 0) THEN
                  CLOSE(ird)
                  CLOSE(iwr)
                  ier_msg(1) = 'Error reading CrystalMaker file'
                  WRITE(ier_msg(2),5000) nline
                  RETURN
               ENDIF
               length = len_str (line) 
               READ(line,1500) latt(1),latt(2),latt(3)
               nline = nline + 1
               IF(iostatus /= 0) THEN
                  CLOSE(ird)
                  CLOSE(iwr)
                  ier_msg(1) = 'Error reading CrystalMaker file'
                  WRITE(ier_msg(2),5000) nline
                  RETURN
               ENDIF
               length = len_str (line) 
               READ(line,1500) latt(1),latt(2),latt(3)
            ELSEIF(str_comp(line(1:4),'List', 4, length, 4)) THEN
               indx1 = INDEX(line,'all')
               indx2 = INDEX(line,'atoms')
               READ(line(indx1+3:indx2-1),*,iostat=iostatus) natoms
               IF(iostatus /= 0) THEN
                  CLOSE(ird)
                  CLOSE(iwr)
                  ier_num    = -3
                  ier_typ    = ER_IO
                  ier_msg(1) = 'Error reading CrystalMaker file'
                  WRITE(ier_msg(2),5000) nline
                  RETURN
               ENDIF
            ELSEIF(str_comp(line(1:4),'Elmt', 4, length, 4)) THEN
               EXIT header
            ENDIF cmd
         ENDIF zero_o
      ENDDO header
      READ (ird, 1000, iostat=iostatus) line
      nline = nline + 1
!
      WRITE (iwr, 2000)         ! Write 'title' line
      WRITE (iwr, 2100)         ! Write 'spcgr P1' line
      WRITE (iwr, 2200) latt    ! Write lattice constants
      WRITE (iwr, 2300)         ! Write 'atoms' line
!
      DO i=1,natoms             ! Loop over all atoms expected in input
         READ (ird, 1000, iostat=iostatus) line
         nline = nline + 1
         IF(iostatus /= 0) THEN
            CLOSE(ird)
            CLOSE(iwr)
            ier_msg(1) = 'Error reading CrystalMaker file'
            WRITE(ier_msg(2),5000) nline
            RETURN
         ENDIF
         WRITE(iwr,3000) line(1:2),line(13:24),line(25:36),line(37:48)
      ENDDO
      CLOSE(ird)
      CLOSE(iwr)
!                                                                       
 1000 FORMAT (a) 
 1500 FORMAT (7x,F10.7,8x,f10.7,9x,f10.7)
 2000 FORMAT ('title ') 
 2100 FORMAT ('spcgr P1') 
 2200 FORMAT ('cell ',5(2x,f9.4,','),2x,f9.4) 
 2300 FORMAT ('atoms') 
 3000 FORMAT (a2,4x,a13,',',a13,',',a13,',  0.100000,   1')
 5000 FORMAT ('Line ',i10)
!                                                                       
!                                                                       
      END SUBROUTINE cmaker2discus                     
!
!*****7*************************************************************************
!
SUBROUTINE rmcprofile2discus (ianz, cpara, lpara, MAXW, ofile, lperiod) 
!-                                                                      
!     converts a RMCProfile "cssr" or "rmc6f" file to DISCUS                   
!+                                                                      
USE build_name_mod
USE precision_mod
USE take_param_mod
USE str_comp_mod
USE support_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER          , INTENT(INOUT)                 :: ianz 
INTEGER          , INTENT(IN)                    :: MAXW 
CHARACTER (LEN=*), DIMENSION(1:MAXW), INTENT(INOUT) :: cpara
INTEGER          , DIMENSION(1:MAXW), INTENT(INOUT) :: lpara
CHARACTER(LEN=*)                    , INTENT(OUT)   :: ofile 
LOGICAL                             , INTENT(IN)    :: lperiod   ! Attempt to rearrange periodically 
!                                                                       
INTEGER, PARAMETER    :: RMC_CSSR  = 0
INTEGER, PARAMETER    :: RMC_RMCF6 = 1
!                                                                       
REAL(KIND=PREC_DP)   , DIMENSION(3) :: werte
!                                                                       
CHARACTER(LEN=PREC_STRING)   :: infile = ' '
INTEGER               :: ird, iwr 
INTEGER               :: style
LOGICAL               :: fileda
!     INTEGER, PARAMETER    :: NOPTIONAL = 1
!     CHARACTER(LEN=   4)       , DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
!     CHARACTER(LEN=PREC_STRING), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
!     INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
!     INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
!     LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para present
!     REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
!     INTEGER, PARAMETER                        :: ncalc = 0 ! Number of values to calculate 
!
!
!     DATA oname  / 'sort'   /
!     DATA loname /  4       /
!     opara  =  (/ 'none' /)   ! Always provide fresh default values
!     lopara =  (/  4     /)
!     owerte =  (/  0.0   /)
!
!     CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
!                       oname, loname, opara, lopara, lpresent, owerte)
!     IF(ier_num/=0) RETURN
!     lperiod = str_comp(opara(1), 'discus', 3, lopara(1), 6)
!
!                                                                       
!     Create input / output file name
!
CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
IF (ier_num.ne.0) THEN 
   RETURN 
ENDIF 
infile = cpara (1)(1:lpara (1) )
INQUIRE(file=infile,exist=fileda)
IF(fileda) THEN
   IF(infile(lpara(1)-4:lpara(1)) == '.cssr' .OR. &
      infile(lpara(1)-4:lpara(1)) == '.CSSR' ) THEN 
      style = RMC_CSSR
      ofile  = cpara (1) (1:lpara (1)-5 ) //'.stru' 
   ELSEIF(infile(lpara(1)-5:lpara(1)) == '.rmc6f' .OR. &
      infile(lpara(1)-5:lpara(1)) == '.RMC6F' ) THEN 
      style = RMC_RMCF6
      ofile  = cpara (1) (1:lpara (1)-6 ) //'.stru' 
   ELSE
      ier_num = -6
      ier_typ = ER_COMM
      RETURN
   ENDIF
ELSE
   infile = cpara (1) (1:lpara (1) ) //'.rmc6f'
   INQUIRE(file=infile,exist=fileda)
   IF(fileda) THEN
      style = RMC_RMCF6
      ofile  = cpara (1) (1:lpara (1) ) //'.stru' 
   ELSE
      infile = cpara (1) (1:lpara (1) ) //'.rmc6f'
      INQUIRE(file=infile,exist=fileda)
      IF(fileda) THEN
         style = RMC_RMCF6
         ofile  = cpara (1) (1:lpara (1) ) //'.stru' 
      ELSE
         infile = cpara (1) (1:lpara (1) ) //'.cssr'
         INQUIRE(file=infile,exist=fileda)
         IF(fileda) THEN
            style = RMC_CSSR
            ofile  = cpara (1) (1:lpara (1) ) //'.stru' 
         ELSE
            infile = cpara (1) (1:lpara (1) ) //'.CSSR'
            INQUIRE(file=infile,exist=fileda)
            IF(fileda) THEN
               style = RMC_CSSR
               ofile  = cpara (1) (1:lpara (1) ) //'.stru' 
            ELSE
               ier_num = -6
               ier_typ = ER_COMM
               RETURN
            ENDIF
         ENDIF
      ENDIF
   ENDIF
ENDIF
!
ird = 34 
iwr = 35 
CALL oeffne (ird, infile, 'old') 
IF (ier_num /= 0) THEN 
   CLOSE(ird)
   RETURN 
ENDIF 
CALL oeffne (iwr, ofile, 'unknown') 
IF(ier_num /= 0) THEN 
   CLOSE(iwr)
   RETURN 
ENDIF 
!                                                                       
IF(style == RMC_CSSR) THEN
   CALL cssr2discus(ird, iwr)
ELSEIF(style == RMC_RMCF6) THEN
   CALL rmcf62discus(ird, iwr, lperiod)
ENDIF
CLOSE(iwr)
CLOSE(ird)
!
END SUBROUTINE rmcprofile2discus 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE cssr2discus(ird, iwr)
!
USE errlist_mod
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: ird
INTEGER, INTENT(IN) :: iwr
!
CHARACTER(LEN= 4)     :: atom   = ' '
CHARACTER(LEN=87)     :: line   = ' '
CHARACTER(LEN=80)     :: title  = ' '
      INTEGER               :: i
INTEGER               :: nline
INTEGER               :: iostatus
INTEGER               :: natoms
REAL   , DIMENSION(6) :: latt! (6) 
REAL   , DIMENSION(3) :: pos ! (6) 
!
      nline     = 1
!
      READ(ird, *    ,IOSTAT=iostatus) latt(1:3)
      IF(iostatus/=0) THEN
         ier_num = -48
         WRITE(ier_msg(1),'(a,i7)') 'Error in line ', nline
         CLOSE(iwr)
         CLOSE(ird)
         RETURN
      ENDIF
      nline     = nline + 1
      READ(ird, *    ,IOSTAT=iostatus) latt(4:6)
      IF(iostatus/=0) THEN
         ier_num = -48
         WRITE(ier_msg(1),'(a,i7)') 'Error in line ', nline
         CLOSE(iwr)
         CLOSE(ird)
         RETURN
      ENDIF
      nline     = nline + 1
      READ(ird, *    ,IOSTAT=iostatus) natoms
      IF(iostatus/=0) THEN
         ier_num = -119
         ier_typ = ER_APPL
         WRITE(ier_msg(1),'(a,i7)') 'Error in line ', nline
         CLOSE(iwr)
         CLOSE(ird)
         RETURN
      ENDIF
      nline     = nline + 1
      READ(ird, '(a)',IOSTAT=iostatus) line
      IF(iostatus/=0) THEN
         ier_num = -46
         ier_typ = ER_APPL
         WRITE(ier_msg(1),'(a,i7)') 'Error in line ', nline
         CLOSE(iwr)
         CLOSE(ird)
         RETURN
      ENDIF
!
      i = INDEX(line,';')
      IF(i > 1) THEN
         title = line(1:i-1)
      ELSE
         title = ' '
      ENDIF
!
      WRITE(iwr, 1000) title
      WRITE(iwr, 1100)
      WRITE(iwr, 1200) latt
      WRITE(iwr, 1300)
!
      atoms: DO i=1,natoms
         nline     = nline + 1
         READ(ird, '(a)',IOSTAT=iostatus) line
         IF(iostatus/=0) THEN
            ier_num = -49
            ier_typ = ER_APPL
            WRITE(ier_msg(1),'(a,i7)') 'Error in line ', nline
            CLOSE(iwr)
            CLOSE(ird)
            RETURN
         ENDIF
         atom = line(8:9)
         READ(line(15:49),*,IOSTAT=iostatus) pos
         IF(iostatus/=0) THEN
            ier_num = -49
            ier_typ = ER_APPL
            WRITE(ier_msg(1),'(a,i7)') 'Error in line ', nline
            CLOSE(iwr)
            CLOSE(ird)
            RETURN
         ENDIF
         WRITE(iwr,2000) atom,pos
      ENDDO atoms
!
1000 FORMAT('title ',a)
1100 FORMAT('spcgr P1')
1200 FORMAT('cell ', 6(2x,F12.6:,', '))
1300 FORMAT('atoms')     
2000 FORMAT(A4,3(2x, F10.6,','),'   0.1,    1')
!
      CLOSE(iwr)
      CLOSE(ird)
!
END SUBROUTINE cssr2discus
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE rmcf62discus(ird, iwr, lperiod)
!
USE blanks_mod
USE errlist_mod
USE get_params_mod
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: ird     ! Input file access number
INTEGER, INTENT(IN) :: iwr     ! Output file access number
LOGICAL, INTENT(IN) :: lperiod ! Attempt to rearrange periodically 
!
CHARACTER(LEN=256)     :: line   = ' '
CHARACTER(LEN=256)     :: title  = ' '
!
CHARACTER(LEN= 4), DIMENSION(:)  , ALLOCATABLE :: r6_at_name
REAL             , DIMENSION(:,:), ALLOCATABLE :: r6_pos
INTEGER          , DIMENSION(:)  , ALLOCATABLE :: r6_site
INTEGER          , DIMENSION(:,:), ALLOCATABLE :: r6_cell
INTEGER, PARAMETER  :: MAXP = 9
INTEGER                               :: ianz
CHARACTER(LEN=1024),DIMENSION(1:MAXP) :: cpara
INTEGER            ,DIMENSION(1:MAXP) :: lpara
INTEGER                               :: length
!
INTEGER                :: i, inumber ! Dummy index
INTEGER                :: iostatus   ! Current line number for error reports
INTEGER                :: nline      ! Current line number for error reports
INTEGER                :: natoms     ! Current line number for error reports
INTEGER                :: nsites     ! Numberr of sites in a unit cell
INTEGER                :: success    !
INTEGER, DIMENSION(3)  :: super      ! Super cell dimensions
REAL                   :: density    ! Number density in Atoms / A^3
REAL   , DIMENSION(6)  :: lattice    ! Unit  cell dimensions For large cell
REAL   , DIMENSION(3,3):: orient     ! Unit  cell dimensions For large cell
!
nline   = 0
natoms  = 0
nsites  = 0
success = 1
header: DO
   nline = nline + 1
   READ(ird,'(a)',iostat=iostatus) line
   IF ( IS_IOSTAT_END(iostatus )) THEN
      WRITE(ier_msg(1),'(a,i8)') 'Error in line ', nline
      ier_num = -6
      ier_typ = ER_IO
      RETURN
   ENDIF
   IF(line(1: 6) == 'Atoms:') EXIT header
   IF(line(1:18) == 'Metadata material:') THEN
      READ(line(19:LEN_TRIM(line)),'(a)') title
      length = LEN_TRIM(title)
      CALL rem_leading_bl(title, length)
   ENDIF
   IF(line(1:16) == 'Number of atoms:') THEN
      READ(line(17:LEN_TRIM(line)),*,IOSTAT=iostatus) natoms
   ENDIF
   IF(line(1:21) == 'Supercell dimensions:') THEN
      READ(line(22:LEN_TRIM(line)),*,IOSTAT=iostatus) super
   ENDIF
   IF(line(1:24) == 'Number density (Ang^-3):') THEN
      READ(line(25:LEN_TRIM(line)),*,IOSTAT=iostatus) density
   ENDIF
   IF(line(1:15) == 'Cell (Ang/deg):') THEN
      READ(line(16:LEN_TRIM(line)),*,IOSTAT=iostatus) lattice
   ENDIF
   IF(line(1:22) == 'Lattice vectors (Ang):') THEN
      nline = nline + 1
      READ(ird,'(a)',iostat=iostatus) line
      READ(line(1:LEN_TRIM(line)),*,IOSTAT=iostatus) orient(1,:)
      nline = nline + 1
      READ(ird,'(a)',iostat=iostatus) line
      READ(line(1:LEN_TRIM(line)),*,IOSTAT=iostatus) orient(2,:)
      nline = nline + 1
      READ(ird,'(a)',iostat=iostatus) line
      READ(line(1:LEN_TRIM(line)),*,IOSTAT=iostatus) orient(3,:)
   ENDIF
ENDDO header
!
!
ALLOCATE(r6_at_name(    1:natoms))
ALLOCATE(r6_pos    (1:3,1:natoms))
ALLOCATE(r6_site   (    1:natoms))
ALLOCATE(r6_cell   (1:3,1:natoms))
r6_at_name(  :) = ' '
r6_pos    (:,:) = 0.0
r6_site   (  :) = 0
r6_cell   (:,:) = 0
atoms:DO i=1,natoms
   nline = nline + 1
   READ(ird,'(a)',iostat=iostatus) line
   length = LEN_TRIM(line)
   CALL get_params_blank (line, ianz, cpara, lpara, MAXP, length)
   IF ( IS_IOSTAT_END(iostatus )) THEN
      WRITE(ier_msg(1),'(a,i8)') 'Error in line ', nline
      ier_num = -6
      ier_typ = ER_IO
      RETURN
   ENDIF
   length = LEN_TRIM(line)
   CALL get_params_blank (line, ianz, cpara, lpara, MAXP, length)
   IF(ier_num/=0) then
      WRITE(ier_msg(1),'(a,i8)') 'Error in line ', nline
      RETURN
   ENDIF
   READ(cpara(1) ,*) inumber
   READ(cpara(2), '(a)') r6_at_name(i)
   READ(cpara(3), *) r6_pos(1,i)
   READ(cpara(4), *) r6_pos(2,i)
   READ(cpara(5), *) r6_pos(3,i)
   READ(cpara(6), *) r6_site(i)
   READ(cpara(7), *) r6_cell(1,i)
   READ(cpara(8), *) r6_cell(2,i)
   READ(cpara(9), *) r6_cell(3,i)
   nsites = MAX(nsites, r6_site(i))
ENDDO atoms
!
IF(lperiod) THEN
   CALL rmc6f_period(natoms, nsites, lattice, super, r6_at_name, r6_pos, r6_site, r6_cell)
ENDIF
WRITE(iwr, 1000) title(1:LEN_TRIM(title))
WRITE(iwr, 1100)
WRITE(iwr, 1200) lattice
IF(lperiod) THEN
   WRITE(iwr, 1250) super, nsites
ELSE
   WRITE(iwr, 1250) 1,1,1, natoms
ENDIF
WRITE(iwr, 1300)
watoms: DO i=1,natoms
   WRITE(iwr, 2000) r6_at_name(i), r6_pos(:,i)
ENDDO watoms
!
DEALLOCATE(r6_at_name)
DEALLOCATE(r6_pos    )
DEALLOCATE(r6_site   )
DEALLOCATE(r6_cell   )
!
1000 FORMAT('title ',a)
1100 FORMAT('spcgr P1')
1200 FORMAT('cell ', 6(2x,F12.6:,', '))
1250 FORMAT('ncell ',3(2x,I6,','),2x,I12)
1300 FORMAT('atoms')     
2000 FORMAT(A4,3(2x, F10.6,','),'   0.1,    1')
!
END SUBROUTINE rmcf62discus
!
!*******************************************************************************
!
SUBROUTINE rmc6f_period(natoms, nsites, lattice, super, r6_at_name, r6_pos, r6_site, r6_cell)
!
INTEGER                          , INTENT(IN)    :: natoms     ! Number of atoms
INTEGER                          , INTENT(IN)    :: nsites     ! Number of sites in teh unit cell
CHARACTER(LEN= 4), DIMENSION(  natoms), INTENT(INOUT) :: r6_at_name
REAL             , DIMENSION(3       ), INTENT(INOUT) :: lattice
INTEGER          , DIMENSION(3       ), INTENT(INOUT) :: super
REAL             , DIMENSION(3,natoms), INTENT(INOUT) :: r6_pos
INTEGER          , DIMENSION(  natoms), INTENT(INOUT) :: r6_site
INTEGER          , DIMENSION(3,natoms), INTENT(INOUT) :: r6_cell
!
REAL, PARAMETER :: EPS = 1.E-6
!
REAL             , DIMENSION(3,2)          :: xyz
REAL             , DIMENSION(3,    nsites) :: shift
REAL             , DIMENSION(3, 2, nsites) :: av_pos
REAL             , DIMENSION(3, 2, nsites) :: si_pos
REAL             , DIMENSION(3,    nsites) :: ave_pos
REAL             , DIMENSION(3,    nsites) :: sig_pos
INTEGER          , DIMENSION(      nsites) :: nav_pos
!
CHARACTER(LEN= 4), DIMENSION(   nsites, super(1), super(2), super(3)) :: dis_atom
REAL             , DIMENSION(3, nsites, super(1), super(2), super(3)) :: dis_pos
!
INTEGER  :: i,j, i1,i2,i3
INTEGER, DIMENSION(3) :: k
REAL   , DIMENSION(3) :: wrap
!
av_pos(:,:,:) = 0.00
si_pos(:,:,:) = 0.00
ave_pos(:,:)  = 0.00
sig_pos(:,:)  = 0.00
nav_pos( :)   = 0
dis_atom(  :,:,:,:) = 'WRNG'
dis_pos (:,:,:,:,:) = 0.0
!
! Check if the number of atoms is an integer multiple of the number of sites
!
IF(ABS(REAL(natoms)/REAL(nsites) - natoms/nsites) > EPS .OR.         &
   super(1)*super(2)*super(3)*nsites /= natoms               ) THEN 
   ier_num = -146
   ier_typ = ER_APPL
   RETURN
ENDIF
!
! Determine average positions, and the shift if atoms a closer to 0,0,0
DO i=1,natoms
   xyz(1,1) = r6_pos(1,i)*super(1)     -REAL( INT(r6_pos(1,i)*super(1)))
   xyz(2,1) = r6_pos(2,i)*super(2)     -REAL( INT(r6_pos(2,i)*super(2)))
   xyz(3,1) = r6_pos(3,i)*super(3)     -REAL( INT(r6_pos(3,i)*super(3)))
   xyz(1,2) = r6_pos(1,i)*super(1)+0.5 -REAL( INT(r6_pos(1,i)*super(1)+0.5))
   xyz(2,2) = r6_pos(2,i)*super(2)+0.5 -REAL( INT(r6_pos(2,i)*super(2)+0.5))
   xyz(3,2) = r6_pos(3,i)*super(3)+0.5 -REAL( INT(r6_pos(3,i)*super(3)+0.5))
!
   av_pos(1,1,r6_site(i)) = av_pos(1,1,r6_site(i)) + xyz(1,1)
   av_pos(2,1,r6_site(i)) = av_pos(2,1,r6_site(i)) + xyz(2,1)
   av_pos(3,1,r6_site(i)) = av_pos(3,1,r6_site(i)) + xyz(3,1)
   av_pos(1,2,r6_site(i)) = av_pos(1,2,r6_site(i)) + xyz(1,2)
   av_pos(2,2,r6_site(i)) = av_pos(2,2,r6_site(i)) + xyz(2,2)
   av_pos(3,2,r6_site(i)) = av_pos(3,2,r6_site(i)) + xyz(3,2)
  nav_pos(  r6_site(i)) =nav_pos(  r6_site(i)) + 1
ENDDO
!
DO i=1,nsites
  IF(nav_pos(i)>0) THEN
      av_pos(:,1,i) = av_pos(:,1,i)/nav_pos(i)
      av_pos(:,2,i) = av_pos(:,2,i)/nav_pos(i)
   ENDIF
ENDDO
!
! Determine a sigma for the two calculation modes in order to decide if atoms are 
! closer to 0,0,0 or closer to 1/2,1/2,1/2
!
DO i=1,natoms
   xyz(1,1) = r6_pos(1,i)*super(1)     -REAL( INT(r6_pos(1,i)*super(1)))
   xyz(2,1) = r6_pos(2,i)*super(2)     -REAL( INT(r6_pos(2,i)*super(2)))
   xyz(3,1) = r6_pos(3,i)*super(3)     -REAL( INT(r6_pos(3,i)*super(3)))
   xyz(1,2) = r6_pos(1,i)*super(1)+0.5 -REAL( INT(r6_pos(1,i)*super(1)+0.5))
   xyz(2,2) = r6_pos(2,i)*super(2)+0.5 -REAL( INT(r6_pos(2,i)*super(2)+0.5))
   xyz(3,2) = r6_pos(3,i)*super(3)+0.5 -REAL( INT(r6_pos(3,i)*super(3)+0.5))
!
   si_pos(1,1,r6_site(i)) = si_pos(1,1,r6_site(i)) + (xyz(1,1)-av_pos(1,1,r6_site(i)))**2
   si_pos(2,1,r6_site(i)) = si_pos(2,1,r6_site(i)) + (xyz(2,1)-av_pos(1,1,r6_site(i)))**2
   si_pos(3,1,r6_site(i)) = si_pos(3,1,r6_site(i)) + (xyz(3,1)-av_pos(1,1,r6_site(i)))**2
   si_pos(1,2,r6_site(i)) = si_pos(1,2,r6_site(i)) + (xyz(1,2)-av_pos(1,2,r6_site(i)))**2
   si_pos(2,2,r6_site(i)) = si_pos(2,2,r6_site(i)) + (xyz(2,2)-av_pos(1,2,r6_site(i)))**2
   si_pos(3,2,r6_site(i)) = si_pos(3,2,r6_site(i)) + (xyz(3,2)-av_pos(1,2,r6_site(i)))**2
ENDDO
!
DO i=1,nsites
   IF(nav_pos(i)>0) THEN
      si_pos(:,1,i) = si_pos(:,1,i)/nav_pos(i)
      si_pos(:,2,i) = si_pos(:,2,i)/nav_pos(i)
   ENDIF
   DO j=1,3
      IF(si_pos(j,1,i) < si_pos(j,2,i)) THEN
         ave_pos(j,i) = av_pos(j,1,i)
         sig_pos(j,i) = si_pos(j,1,i)
         shift(j,i)   = 0.0
      ELSE
         ave_pos(j,i) = av_pos(j,2,i) -0.5
         sig_pos(j,i) = si_pos(j,2,i)
         shift(j,i)   = 0.5
      ENDIF
   ENDDO
!write(*,'(a,i2,2x,6(f6.3,2x),i6, 3(f6.3:,2x))') 'Average ',i, ave_pos(:,i), sig_pos(:,i),nav_pos(i),shift(:,i)
ENDDO
! Copy atoms into DISCUS sequence
!
DO i=1,natoms
   DO j=1,3
      k(j) =  INT(r6_pos(j,i)*super(1) + shift(j,r6_site(i))) + 1
      wrap(j) = 0.0
      IF(k(j)<1) THEN
         wrap(j) = super(j)
         k(j)    = k(j) + super(j)
      ELSEIF(k(j)>super(j)) THEN
         wrap(j) = -super(j)
         k(j)    = k(j) -super(j)
      ENDIF
   ENDDO
!
   dis_atom(  r6_site(i),k(1),k(2),k(3)) = r6_at_name(i)
   dis_pos (1,r6_site(i),k(1),k(2),k(3)) = r6_pos(1,  i)*super(1) + wrap(1)
   dis_pos (2,r6_site(i),k(1),k(2),k(3)) = r6_pos(2,  i)*super(2) + wrap(2)
   dis_pos (3,r6_site(i),k(1),k(2),k(3)) = r6_pos(3,  i)*super(3) + wrap(3)
ENDDO
!
! Copy back into linear array 
i=0
DO i3=1,super(3)
   DO i2=1,super(2)
      DO i1=1,super(1)
         DO  j=1,nsites
            i = i + 1
            r6_at_name(i) = dis_atom(  j,i1,i2,i3)
            r6_pos(1,  i) = dis_pos (1,j,i1,i2,i3)
            r6_pos(2,  i) = dis_pos (2,j,i1,i2,i3)
            r6_pos(3,  i) = dis_pos (3,j,i1,i2,i3)
         ENDDO
      ENDDO
   ENDDO
ENDDO
lattice(1:3) = lattice(1:3)/REAL(super(1:3))
!
END SUBROUTINE rmc6f_period
!
!
      SUBROUTINE cif2discus (ianz, cpara, lpara, MAXW, ofile) 
!-                                                                      
!     converts a CIF file to DISCUS                   
!+                                                                      
!                                                                       
!      USE tensors_mod
      USE build_name_mod
      USE wink_mod
      USE ber_params_mod
      USE blanks_mod
      USE get_params_mod
use matrix_mod
USE lib_length
USE precision_mod
      USE string_convert_mod
USE support_mod
!
      IMPLICIT none 
!                                                                       
      INTEGER          , INTENT(INOUT)                 :: ianz 
      INTEGER          , INTENT(IN)                    :: MAXW 
      CHARACTER (LEN=*), DIMENSION(1:MAXW), INTENT(INOUT) :: cpara
      INTEGER          , DIMENSION(1:MAXW), INTENT(INOUT) :: lpara
      CHARACTER(LEN=*)                 , INTENT(OUT)   :: ofile 
!                                                                       
      REAL(KIND=PREC_DP), PARAMETER :: eightpi2 = 8.D0*3.1415926535897932384626433832795028841971693993751D0**2
      REAL, PARAMETER :: EPS = 0.00001
!                                                                       
      REAL(KIND=PREC_DP)   , DIMENSION(3) :: werte
!                                                                       
      CHARACTER(LEN= 1)     :: bravais= ' '
      CHARACTER(LEN=80)     :: title  = ' '
      CHARACTER(LEN=80)     :: newtitle  = ' '
      CHARACTER(LEN=80)     :: spcgr  = ' '
      CHARACTER(LEN=80)     :: aniso_label  = ' '
      CHARACTER(LEN=80)     :: aniso_symb   = ' '
      CHARACTER(LEN=PREC_STRING)   :: infile = ' '
      CHARACTER(LEN=PREC_STRING)   :: wfile  = ' '
      CHARACTER(LEN=PREC_STRING)                              :: line
      CHARACTER(LEN=PREC_STRING)                              :: line_cap
      CHARACTER(LEN=PREC_STRING), DIMENSION(:), ALLOCATABLE   :: rawline
      CHARACTER(LEN=PREC_STRING), DIMENSION(:), ALLOCATABLE   :: ccpara
      CHARACTER(LEN=PREC_STRING), DIMENSION(3)                :: cspara
      INTEGER            , DIMENSION(:), ALLOCATABLE   :: llpara
      INTEGER            , DIMENSION(3)                :: lspara
      REAL(KIND=PREC_DP) , DIMENSION(3)                :: wwerte
      INTEGER               :: MAXLINES 
      INTEGER               :: ird, iwr 
      INTEGER               :: iianz = 3
      INTEGER               :: i, j, k
      INTEGER               :: iblank
      INTEGER               :: iostatus
      LOGICAL, DIMENSION(7) :: header_done = .false.
      INTEGER               :: line_no, line_sig, data_no
      INTEGER               :: length, length_cap, length_b
      INTEGER               :: is_cell
      INTEGER               :: is_loop
      INTEGER               :: is_spcgr
      INTEGER               :: is_spcgr_no
      INTEGER               :: is_symm
      INTEGER               :: is_atom
      INTEGER               :: is_anis
      INTEGER               :: is_paren
      INTEGER               :: j_atom  = 0
      INTEGER               :: j_anis  = 0
      INTEGER               :: j_symb  = 0
      INTEGER               :: j_label = 0
      INTEGER               :: j_uiso  = 0
      INTEGER               :: j_biso  = 0
      INTEGER               :: j_occ   = 0
      INTEGER               :: j_x     = 0
      INTEGER               :: j_y     = 0
      INTEGER               :: j_z     = 0
      INTEGER               :: j_aniso_symb  = 0
      INTEGER               :: j_aniso_label = 0
      INTEGER               :: j_aniso_11    = 0
      INTEGER               :: j_aniso_22    = 0
      INTEGER               :: j_aniso_33    = 0
      INTEGER               :: j_aniso_12    = 0
      INTEGER               :: j_aniso_13    = 0
      INTEGER               :: j_aniso_23    = 0
      INTEGER               :: j_aniso_B11   = 0
      INTEGER               :: j_aniso_B22   = 0
      INTEGER               :: j_aniso_B33   = 0
      INTEGER               :: j_aniso_B12   = 0
      INTEGER               :: j_aniso_B13   = 0
      INTEGER               :: j_aniso_B23   = 0
      INTEGER               :: symm_1, symm_2, symm_n
      INTEGER               :: ix,iy,iz
      INTEGER               :: nentries
      INTEGER               :: spcgr_no
      INTEGER               :: iquote1
      INTEGER               :: iquote2
      INTEGER               :: spcgr_l
      INTEGER               :: nline
      INTEGER               :: nblank
      LOGICAL               :: in_section
      LOGICAL               :: l_space_group
      LOGICAL               :: l_numbers
      INTEGER               :: data_i
      REAL   , DIMENSION(6) :: latt! (6) 
      REAL   , DIMENSION(3) :: pos ! (6) 
      REAL   , DIMENSION(3) :: rlatt    ! (6) 
      REAL   , DIMENSION(3,3) :: uij ! (6) 
      REAL   , DIMENSION(3,3) :: bij ! (6) 
      REAL(kind=PREC_DP)   , DIMENSION(3,3) :: gten ! (6) 
      REAL(kind=PREC_DP)   , DIMENSION(3,3) :: rten ! (6) 
      REAL   , DIMENSION(4,4) :: symm_mat ! (6) 
      REAL                  :: uiso
      REAL                  :: biso
      REAL                  :: occ
!
      TYPE :: atom_list
         CHARACTER (LEN=80) :: label  
         CHARACTER (LEN=80) :: symbol  
         CHARACTER (LEN=4)  :: at_name
         REAL,DIMENSION(3)  :: at_pos
         REAL,DIMENSION(6)  :: at_uanis
         REAL               :: at_bvalue
         REAL               :: at_occ
         TYPE(atom_list), POINTER   :: next
      END TYPE atom_list
!
      TYPE(atom_list), POINTER :: head
      TYPE(atom_list), POINTER :: tail
      TYPE(atom_list), POINTER :: temp
!
!
      is_loop = 0
      symm_n  = 0
      symm_1  = 0
!                                                                       
!     Create input / output file name
!
      CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
      IF (ier_num.ne.0) THEN 
         RETURN 
      ENDIF 
      infile = cpara (1) 
      i = index (infile, '.',.true.)                  ! find last occurence of '.'
      IF (i.eq.0) THEN 
         infile = cpara (1) (1:lpara (1) ) //'.cif' 
         ofile  = cpara (1) (1:lpara (1) ) //'.stru' 
      ELSE 
         IF(    cpara(1)(lpara(1)-3:lpara(1)) == '.cif') THEN
            ofile  = cpara (1) (1:lpara(1)-3) //'stru' 
         ELSEIF(cpara(1)(lpara(1)-3:lpara(1)) == '.CIF') THEN
            ofile  = cpara (1) (1:lpara(1)-3) //'stru' 
         ELSE
            ofile  = cpara (1) (1:i) //'stru'
         ENDIF
      ENDIF 
      ird = 34 
      iwr = 35 
      CALL oeffne (ird, infile, 'old') 
      IF (ier_num.ne.0) THEN 
         RETURN 
      ENDIF 
!
      NULLIFY(head)
      NULLIFY(tail)
      NULLIFY(temp)
!
! As we do not know the length of the input file, lets read it once
!
      line_sig= 0
      line_no = 0
      data_no = 0      ! Counter for individual "data_" sections
countline: DO
         READ(ird, '(a)', IOSTAT=iostatus) line
         IF ( IS_IOSTAT_END(iostatus )) EXIT countline
         line_no = line_no + 1
         length  = len_str(line)
         IF(length > 0 ) THEN
            line_sig = line_no
            CALL rem_leading_bl(line,length)
            CALL do_cap(line)
            IF(line(1:5) == 'DATA_') data_no = data_no + 1
         ENDIF
      ENDDO countline
      MAXLINES = line_sig
      ALLOCATE(rawline(1:MAXLINES))
      rawline = ' '
      line_no = 0
      REWIND(ird)
getline: DO
         READ(ird, '(a)', IOSTAT=iostatus) rawline(line_no+1)
         IF ( IS_IOSTAT_END(iostatus )) EXIT getline
         line_no = line_no + 1
         length  = len_str(rawline(line_no))
         CALL rem_leading_bl(rawline(line_no),length)
         CALL tab2blank(rawline(line_no), length)
         IF(line_no == line_sig) EXIT getline
      ENDDO getline
      CLOSE(ird)
!
      nline     = 0
      in_section = .false.
!
!   Loop over all observed "data_" sections, write to separate files
      data_i = 0
data_entries: DO WHILE(data_i < data_no)
!
!   Run an loop over all input lines
!
main: DO 
         nline = nline + 1
         IF(nline==line_no) EXIT main   ! End of input
         line = rawline(nline)
         length = len_str(line)
         length_cap = length
         line_cap = line
         CALL do_cap(line_cap)
         IF(length   == 0 ) CYCLE main
         IF(line(1:1)=='#') CYCLE main
!
!  Data statement
!
         IF(INDEX(line_cap(1:5),'DATA_')/=0) THEN        ! Found a "data_" line
                  IF(length > 5) THEN
                     newtitle = line(6:length)
                  ENDIF
            IF(in_section) THEN           ! invalid until we find first "data_" line
               IF(data_i < data_no) THEN  ! For all but last branch out to write previous section
                  EXIT main               ! End of previous "data_" section, write file
               ELSE
                  CONTINUE                ! Never reached, as last section has its "data_" read previously
               ENDIF
            ELSE                          ! At the first "data_" line
               title = newtitle           ! immediately save title for write
               in_section = .TRUE.        ! We are now in a "data_" section
            ENDIF
         ENDIF
!
!  Loop statement
!
         IF(INDEX(line_cap,'LOOP_')/=0) THEN
            is_loop = nline                     ! Store line number of loop start
         ENDIF
!
!  Space group name                    ! Find '_H-M' or '_H-M_alt' or similar
!
         is_spcgr = INDEX(line,'_space_group')
         IF(is_spcgr/=0) THEN                   ! Got a symmetry info
!
!           IF(INDEX(line,'_symmetry_space_group_name_H-M')/=0) THEN
            is_spcgr = INDEX(line,'_H-M')
            IF(is_spcgr/=0) THEN                       ! Found H-M symbol
               
               iblank = INDEX(line(is_spcgr:length), ' ')    !Find blank after H-M
                  
               iquote1 = INDEX(line(is_spcgr+iblank:length),'''')
               IF(iquote1>0) THEN
                  iquote2 = INDEX(line(is_spcgr+iblank+iquote1:length),'''')
                  iquote1 = is_spcgr+iblank+iquote1
                  iquote2 =             iquote1+iquote2 - 2
               ELSE
                  iquote1 = INDEX(line(is_spcgr+iblank:length),'"')
                  iquote2 = INDEX(line(is_spcgr+iblank+iquote1:length),'"')
                  iquote1 = is_spcgr+iblank+iquote1
                  iquote2 =             iquote1+iquote2 - 2
               ENDIF
               IF(iquote2> iquote1 .and. iquote1>0 .and. iquote2>0 ) THEN
                  spcgr   = line(iquote1:iquote2)
                  spcgr_l = iquote2 - iquote1 + 1
               ELSE                     ! Space group is not enclosed in quotation marks
                  spcgr = line(is_spcgr+iblank:length)
                  spcgr_l = length - (is_spcgr+iblank) + 1
               ENDIF
               CALL rem_bl(spcgr,spcgr_l)
               bravais = spcgr(1:1)
               CALL do_low(spcgr)        ! Make lower case
               CALL do_cap(bravais)      ! Upper case lattice type
               spcgr(1:1) = bravais
               IF(spcgr(3:3)=='3' .AND. spcgr(2:2)/='-'  &
                                  .AND. spcgr(2:2)/='6') THEN
                  spcgr = spcgr(1:2)//'-'//spcgr(3:spcgr_l)
                  spcgr_l = spcgr_l + 1
               ENDIF
               IF(spcgr(spcgr_l-1:spcgr_l-1)==':') THEN
                  spcgr(spcgr_l-1:spcgr_l-1) =','
               ENDIF
               header_done(1) = .true.
            ENDIF
         ENDIF
!
!  space group number
!
         is_spcgr_no = INDEX(line,'_space_group_IT_number')
         IF(is_spcgr_no/=0) THEN
               READ(line(is_spcgr_no+23:length),*,IOSTAT=iostatus) spcgr_no
               header_done(1) = .true.
         ENDIF
!
!  Symmetry operators 
!
         is_symm  = MAX(INDEX(line,'_symmetry_equiv_pos_as_xyz'),       &
                        INDEX(line,'_space_group_symop_operation_xyz'))
         IF(is_symm/=0) THEN                   ! Got a symmetry info
            symm_1 = nline + 1
            symm_2 = nline
            i = 1
            count_symm: DO
               IF(rawline(nline+i)== ' ') EXIT count_symm
               IF(rawline(nline+i)(1:4)== 'loop') EXIT count_symm
               IF(rawline(nline+i)(1:1)== '_'  ) EXIT count_symm
               IF(rawline(nline+i)(1:1)== ';'  ) EXIT count_symm
               i = i+1
            ENDDO count_symm
            symm_2 = nline + i - 1
            symm_n = symm_2-symm_1+1
         ENDIF
!
!
!  Unit cell dimensions
!
         is_cell = INDEX(line,'_cell_')         ! got a cell info
         IF(is_cell/=0) THEN
            is_paren = INDEX(line,'(')
            IF(is_paren > 0 ) THEN
               length = is_paren-1
            ENDIF
            IF(INDEX(line,'_cell_length_a')/=0) THEN
               READ(line(is_cell+14:length),*,IOSTAT=iostatus) latt(1)
               header_done(2) = .true.
            ELSEIF(INDEX(line,'_cell_length_b')/=0) THEN
               READ(line(is_cell+14:length),*,IOSTAT=iostatus) latt(2)
               header_done(3) = .true.
            ELSEIF(INDEX(line,'_cell_length_c')/=0) THEN
               READ(line(is_cell+14:length),*,IOSTAT=iostatus) latt(3)
               header_done(4) = .true.
            ELSEIF(INDEX(line,'_cell_angle_alpha')/=0) THEN
               READ(line(is_cell+17:length),*,IOSTAT=iostatus) latt(4)
               header_done(5) = .true.
            ELSEIF(INDEX(line,'_cell_angle_beta')/=0) THEN
               READ(line(is_cell+16:length),*,IOSTAT=iostatus) latt(5)
               header_done(6) = .true.
            ELSEIF(INDEX(line,'_cell_angle_gamma')/=0) THEN
               READ(line(is_cell+17:length),*,IOSTAT=iostatus) latt(6)
               header_done(7) = .true.
            ENDIF
         ENDIF
!
!  atom coordinates
!
         is_atom = INDEX(line,'_atom_site_fract_x')
         IF(is_atom /= 0) THEN               ! found the loop with atom coordinates
            j_atom = is_loop                 ! start in line after 'loop_'
            nentries = 0
            j_label = 0
            j_symb  = 0
            j_uiso  = 0
            j_biso  = 0
            j_occ   = 0
            j_x     = 0
            j_y     = 0
            j_z     = 0
analyze_atom: DO
               j_atom = j_atom + 1
               line = rawline(j_atom)
               length = len_str(line)
               IF(line(1:1)=='#' .or. line == ' ') CYCLE analyze_atom
               IF(line(1:10)/='_atom_site' .or. j_atom>line_no) THEN
                  IF(j_atom < nline) THEN      ! wrong line prior to '_atom_site_frac_x'
                     EXIT main
                  ENDIF
                  nline = j_atom             ! We are now in line j_atom
                  EXIT analyze_atom
               ENDIF
               nentries = nentries + 1
               IF(line(1:16)=='_atom_site_label')           j_label = nentries
               IF(line(1:22)=='_atom_site_type_symbol')     j_symb  = nentries
               IF(line(1:25)=='_atom_site_U_iso_or_equiv')  j_uiso  = nentries
               IF(line(1:25)=='_atom_site_B_iso_or_equiv')  j_biso  = nentries
               IF(line(1:20)=='_atom_site_occupancy')       j_occ   = nentries
               IF(line(1:18)=='_atom_site_fract_x')         j_x     = nentries
               IF(line(1:18)=='_atom_site_fract_y')         j_y     = nentries
               IF(line(1:18)=='_atom_site_fract_z')         j_z     = nentries
            ENDDO analyze_atom
            IF(.NOT. ALLOCATED(CCPARA)) ALLOCATE(ccpara(nentries))
            IF(.NOT. ALLOCATED(LLPARA)) ALLOCATE(llpara(nentries))
            ccpara = ' '
            nblank = 0
atoms:      DO                                 ! Get all atoms information
               IF(line(1:4)=='loop') THEN 
                  is_loop = nline
                  EXIT atoms
               ENDIF
               IF(line(1:1)/='#' .and. line /= ' ') THEN
                  CALL get_params_blank(line,ianz, ccpara,llpara, nentries, length)
!
!   If there are a different number of parameters, the line does not appear to be 
!   another atom line
!
                  IF(nentries/=ianz) THEN         ! no more atom lines
                     nline = j_atom
                     IF(line(1:4)=='loop') THEN 
                        is_loop = nline
                     ENDIF
                     EXIT atoms
                  ENDIF
                  DO i=1,ianz                     ! Replace all '.' by '0.0'
                     IF(ccpara(i) == '.') THEN
                        ccpara(i) = '0.0'
                        llpara(i) = 3
                     ENDIF
                  ENDDO
                  IF(INDEX(ccpara(j_x),'(')>0) llpara(j_x) = INDEX(ccpara(j_x),'(') - 1
                  IF(INDEX(ccpara(j_y),'(')>0) llpara(j_y) = INDEX(ccpara(j_y),'(') - 1
                  IF(INDEX(ccpara(j_z),'(')>0) llpara(j_z) = INDEX(ccpara(j_z),'(') - 1
                  READ(ccpara(j_x)(1:llpara(j_x)),*) pos(1)
                  READ(ccpara(j_y)(1:llpara(j_y)),*) pos(2)
                  READ(ccpara(j_z)(1:llpara(j_z)),*) pos(3)
                  uiso = 0.0
                  biso = 0.0
                  occ  = 1.0
                  IF(j_uiso > 0) THEN
                     IF(ccpara(j_uiso)(1:1)=='?') THEN
                       biso = 0.000
                     ELSE
                        IF(ccpara(j_uiso)(1:llpara(j_uiso))=='.') THEN
                           biso = 0.0
                        ELSE
                           IF(INDEX(ccpara(j_uiso),'(')>0) llpara(j_uiso) = INDEX(ccpara(j_uiso),'(') - 1
                           READ(ccpara(j_uiso)(1:llpara(j_uiso)),*) uiso
                           biso = uiso * eightpi2
                        ENDIF
                     ENDIF
                  ELSEIF(j_biso > 0) THEN
                     IF(ccpara(j_biso)(1:1)=='?') THEN
                       biso = 0.000
                     ELSE
                        IF(ccpara(j_biso)(1:llpara(j_biso))=='.') THEN
                           biso = 0.0
                        ELSE
                           IF(INDEX(ccpara(j_biso),'(')>0) llpara(j_biso) = INDEX(ccpara(j_biso),'(') - 1
                           READ(ccpara(j_biso)(1:llpara(j_biso)),*) biso
                        ENDIF
                     ENDIF
                  ENDIF
                  IF(j_occ  > 0) THEN
                     IF(ccpara(j_occ )(1:1)=='?') THEN
                       occ  = 1.000
                     ELSE
                        IF(INDEX(ccpara(j_occ ),'(')>0) llpara(j_occ ) = INDEX(ccpara(j_occ ),'(') - 1
                        READ(ccpara(j_occ )(1:llpara(j_occ )),*) occ 
                     ENDIF
                  ENDIF
                  ALLOCATE(TEMP)
                  TEMP%at_name   = ' '
                  TEMP%at_bvalue = 0.0
                  TEMP%at_occ    = 0.0
                  TEMP%at_pos    = 0.0
                  IF(j_symb > 0) THEN  ! I prefer the atom symbol to its label
                     TEMP%label     = ccpara(j_label)(1:      llpara(j_label) )
                     TEMP%symbol    = ccpara(j_symb )(1:      llpara(j_symb ) )
                     TEMP%at_name   = ccpara(j_symb )(1:MIN(4,llpara(j_symb )))
                     CALL test_atom_name(.TRUE.,TEMP%at_name)
                  ELSE
                     TEMP%label     = ccpara(j_label)(1:      llpara(j_label) )
                     TEMP%symbol    = ' '
                     TEMP%at_name   = ccpara(j_label)(1:MIN(4,llpara(j_label)))
                     CALL test_atom_name(.FALSE.,TEMP%at_name)
                  ENDIF
                  TEMP%at_pos(1) = pos(1)
                  TEMP%at_pos(2) = pos(2)
                  TEMP%at_pos(3) = pos(3)
                  TEMP%at_bvalue = biso
                  TEMP%at_occ    = occ
                  TEMP%at_uanis  = 0.0
                  NULLIFY(temp%next)
!
                  IF(ASSOCIATED(TAIL)) THEN
                     TAIL%NEXT => TEMP
                     TAIL      => TAIL%NEXT
                  ELSE
                     TAIL      => TEMP
                     HEAD      => TEMP
                  ENDIF
                  j_atom = j_atom + 1
               ELSE    ! Comment or empty line
                  nblank = nblank + 1
               ENDIF   ! end no comment
!
               IF(j_atom+nblank > line_no) THEN
                  nline = j_atom + nblank            ! We are now in line j_atom
                  IF(ALLOCATED(CCPARA)) DEALLOCATE(ccpara)
                  IF(ALLOCATED(LLPARA)) DEALLOCATE(llpara)
                  EXIT main
               ENDIF
               line   = rawline(j_atom + nblank)
               length = len_str(line)
               nline = j_atom + nblank            ! We are now in line j_atom
            ENDDO atoms
            IF(ALLOCATED(CCPARA)) DEALLOCATE(ccpara)
            IF(ALLOCATED(LLPARA)) DEALLOCATE(llpara)
         ENDIF
!
!  anisotropic displacement parameters
!
         is_anis = INDEX(line,'_atom_site_aniso')
         IF(is_anis /= 0) THEN               ! found the loop with aniso ADP's
            j_anis = is_loop                 ! start in line after 'loop_'
         ENDIF
      ENDDO main
!
!  The main atom list did not contain isotropic U/B values, 
!  obtain equivalent values from the anisotropic ADP's
!
      IF(j_uiso == 0 .AND. j_anis > 0) THEN  ! Found anisotropic ADP's, and no ISO
         nline    = j_anis
         nentries = 0
analyze_anis: DO
            nline = nline + 1
            line = rawline(nline)
            length = len_str(line)
            IF(line(1:1)/='#')  THEN
               IF(line(1:16)/='_atom_site_aniso' .or. nline>line_no) THEN
                  EXIT analyze_anis
               ENDIF
               nentries = nentries + 1
               IF(line(1:22)=='_atom_site_aniso_label')           j_aniso_label = nentries
               IF(line(1:28)=='_atom_site_aniso_type_symbol')     j_aniso_symb  = nentries
               IF(line(1:21)=='_atom_site_aniso_U_11')            j_aniso_11    = nentries
               IF(line(1:21)=='_atom_site_aniso_U_22')            j_aniso_22    = nentries
               IF(line(1:21)=='_atom_site_aniso_U_33')            j_aniso_33    = nentries
               IF(line(1:21)=='_atom_site_aniso_U_12')            j_aniso_12    = nentries
               IF(line(1:21)=='_atom_site_aniso_U_13')            j_aniso_13    = nentries
               IF(line(1:21)=='_atom_site_aniso_U_23')            j_aniso_23    = nentries
               IF(line(1:21)=='_atom_site_aniso_B_11')            j_aniso_B11   = nentries
               IF(line(1:21)=='_atom_site_aniso_B_22')            j_aniso_B22   = nentries
               IF(line(1:21)=='_atom_site_aniso_B_33')            j_aniso_B33   = nentries
               IF(line(1:21)=='_atom_site_aniso_B_12')            j_aniso_B12   = nentries
               IF(line(1:21)=='_atom_site_aniso_B_13')            j_aniso_B13   = nentries
               IF(line(1:21)=='_atom_site_aniso_B_23')            j_aniso_B23   = nentries
            ENDIF
         ENDDO analyze_anis
!
!           Build metric tensors and get recipr. lattice params
!
         gten(1,1) = latt(1)**2
         gten(2,2) = latt(2)**2
         gten(3,3) = latt(3)**2
         gten(1,2) = latt(1)*latt(2)*cos(REAL(rad)*latt(6))
         gten(1,3) = latt(1)*latt(3)*cos(REAL(rad)*latt(5))
         gten(2,3) = latt(2)*latt(3)*cos(REAL(rad)*latt(4))
         gten(2,1) = gten(1,2)
         gten(3,1) = gten(1,3)
         gten(2,3) = gten(3,2)
!        CALL invmat(rten,gten)
         CALL matinv(gten,rten)
         rlatt(1) = SQRT(rten(1,1))
         rlatt(2) = SQRT(rten(2,2))
         rlatt(3) = SQRT(rten(3,3))
         IF(ALLOCATED(CCPARA)) DEALLOCATE(ccpara)
         IF(ALLOCATED(LLPARA)) DEALLOCATE(llpara)
         zero_entries:IF(nentries>0) THEN
         IF(.NOT. ALLOCATED(CCPARA)) ALLOCATE(ccpara(nentries))
         IF(.NOT. ALLOCATED(LLPARA)) ALLOCATE(llpara(nentries))
         ccpara = ' '
anis:    DO                                 ! Get all anisotropic information
            IF(line(1:1)/='#')  THEN
            CALL get_params_blank(line,ianz, ccpara,llpara, nentries, length)
!
!   If there are a different number of parameters, the line does not appear to be 
!   another atom line
!
            IF(nentries/=ianz) THEN         ! no more aniso lines
               EXIT anis
            ENDIF
            aniso_label = ' '
            aniso_symb  = ' '
            IF(j_aniso_label > 0) THEN
               aniso_label = ccpara(j_aniso_label)(1:llpara(j_aniso_label))
            ENDIF
            IF(j_aniso_symb  > 0) THEN
               aniso_symb  = ccpara(j_aniso_symb )(1:llpara(j_aniso_symb ))
            ENDIF
            IF(j_aniso_11 > 0 ) THEN
               IF(INDEX(ccpara(j_aniso_11),'(')>0) llpara(j_aniso_11) = INDEX(ccpara(j_aniso_11),'(') - 1
               IF(INDEX(ccpara(j_aniso_22),'(')>0) llpara(j_aniso_22) = INDEX(ccpara(j_aniso_22),'(') - 1
               IF(INDEX(ccpara(j_aniso_33),'(')>0) llpara(j_aniso_33) = INDEX(ccpara(j_aniso_33),'(') - 1
               IF(INDEX(ccpara(j_aniso_12),'(')>0) llpara(j_aniso_12) = INDEX(ccpara(j_aniso_12),'(') - 1
               IF(INDEX(ccpara(j_aniso_13),'(')>0) llpara(j_aniso_13) = INDEX(ccpara(j_aniso_13),'(') - 1
               IF(INDEX(ccpara(j_aniso_23),'(')>0) llpara(j_aniso_23) = INDEX(ccpara(j_aniso_23),'(') - 1
               READ(ccpara(j_aniso_11)(1:llpara(j_aniso_11)),*) uij(1,1)
               READ(ccpara(j_aniso_22)(1:llpara(j_aniso_22)),*) uij(2,2)
               READ(ccpara(j_aniso_33)(1:llpara(j_aniso_33)),*) uij(3,3)
               READ(ccpara(j_aniso_12)(1:llpara(j_aniso_12)),*) uij(1,2)
               READ(ccpara(j_aniso_13)(1:llpara(j_aniso_13)),*) uij(1,3)
               READ(ccpara(j_aniso_23)(1:llpara(j_aniso_23)),*) uij(2,3)
               uij(2,1) = uij(1,2)
               uij(3,1) = uij(1,3)
               uij(3,2) = uij(2,3)
               uiso = 0.0
               DO i=1,3
                  DO j=1,3
                     uiso = uiso + uij(i,j)*latt(i)*latt(j)*rlatt(i)*rlatt(j)
                  ENDDO
               ENDDO
               uiso = uiso / 3.
               biso = uiso * eightpi2
            ELSEIF(j_aniso_B11 > 0 ) THEN
               IF(INDEX(ccpara(j_aniso_B11),'(')>0) llpara(j_aniso_B11) = INDEX(ccpara(j_aniso_B11),'(') - 1
               IF(INDEX(ccpara(j_aniso_B22),'(')>0) llpara(j_aniso_B22) = INDEX(ccpara(j_aniso_B22),'(') - 1
               IF(INDEX(ccpara(j_aniso_B33),'(')>0) llpara(j_aniso_B33) = INDEX(ccpara(j_aniso_B33),'(') - 1
               IF(INDEX(ccpara(j_aniso_B12),'(')>0) llpara(j_aniso_B12) = INDEX(ccpara(j_aniso_B12),'(') - 1
               IF(INDEX(ccpara(j_aniso_B13),'(')>0) llpara(j_aniso_B13) = INDEX(ccpara(j_aniso_B13),'(') - 1
               IF(INDEX(ccpara(j_aniso_B23),'(')>0) llpara(j_aniso_B23) = INDEX(ccpara(j_aniso_B23),'(') - 1
               READ(ccpara(j_aniso_B11)(1:llpara(j_aniso_B11)),*) bij(1,1)
               READ(ccpara(j_aniso_B22)(1:llpara(j_aniso_B22)),*) bij(2,2)
               READ(ccpara(j_aniso_B33)(1:llpara(j_aniso_B33)),*) bij(3,3)
               READ(ccpara(j_aniso_B12)(1:llpara(j_aniso_B12)),*) bij(1,2)
               READ(ccpara(j_aniso_B13)(1:llpara(j_aniso_B13)),*) bij(1,3)
               READ(ccpara(j_aniso_B23)(1:llpara(j_aniso_B23)),*) bij(2,3)
               bij(2,1) = bij(1,2)
               bij(3,1) = bij(1,3)
               bij(3,2) = bij(2,3)
               biso = 0.0
               DO i=1,3
                  DO j=1,3
                     biso = biso + bij(i,j)*latt(i)*latt(j)*rlatt(i)*rlatt(j)
                  ENDDO
               ENDDO
               biso = biso / 3.
            ENDIF
            TEMP => HEAD
find:       DO WHILE (ASSOCIATED(TEMP))
               IF(j_label > 0 .AND. j_aniso_label > 0) THEN
                  IF(TEMP%label == aniso_label) THEN
                     TEMP%at_bvalue = biso
                  ENDIF
               ELSEIF(j_symb > 0 .AND. j_aniso_symb > 0) THEN
                  IF(TEMP%symbol == aniso_symb) THEN
                     TEMP%at_bvalue = biso
                  ENDIF
               ENDIF
               TEMP => TEMP%next
            ENDDO find
            ENDIF   ! no comment
!
            nline = nline + 1
            IF(nline>line_no) THEN
               EXIT anis
            ENDIF
            line   = rawline(nline)
            length = len_str(line)
         ENDDO anis   
         ENDIF zero_entries
         IF(ALLOCATED(CCPARA)) DEALLOCATE(ccpara)
         IF(ALLOCATED(LLPARA)) DEALLOCATE(llpara)
      ENDIF
!
!  Finally, write the structure to file
!
      IF(data_i > 0) THEN                    ! For all but first file append a number
         WRITE(line(1:6),'(I6.6)')  data_i
         wfile= ofile(1:len_str(ofile))//LINE(1:6)
      ELSE                                   ! This is the first file
         wfile = ofile
      ENDIF
      CALL oeffne (iwr, wfile, 'unknown') 
      IF (ier_num.ne.0) THEN       ! Error opening file, clear memory structure
         TAIL => HEAD
         TEMP => HEAD
         DO WHILE (ASSOCIATED(TAIL))
            TAIL => TAIL%next
            DEALLOCATE(TEMP)       ! Clean up the memory structure
            TEMP => TAIL
         ENDDO
         NULLIFY(HEAD)
         NULLIFY(TEMP)
         NULLIFY(TAIL)
         RETURN 
      ENDIF 
      WRITE(iwr, 1000) title(1:len_str(title))
      l_space_group = .FALSE.
      IF(spcgr /= ' ') THEN
         length = LEN_TRIM(spcgr)
         IF(spcgr(1:1) == '?' .OR. spcgr(2:2) == '?') THEN  !'HM is a '?'
            IF(symm_n>0) THEN 
               spcgr = 'P1'
            ELSE                     !, flag error but finish writing
               ier_num = -126
               ier_typ = ER_APPL
            ENDIF
         ELSE
            IF(length > 1) THEN
               l_space_group = spcgr_test(spcgr ) ! Test for known space group
               IF(.NOT. l_space_group) THEN
                  IF(spcgr(2:2)=='1' .AND. spcgr(length:length)=='1') THEN
                     spcgr = spcgr(1:1) // spcgr(3:length-1)
                     l_space_group = spcgr_test(spcgr ) ! Test for known space group
                  ELSE
                     IF(ABS(latt(4)-90.0)<EPS .AND. ABS(latt(6)-90.0)<EPS.AND.  &
                        ABS(latt(5)-90.0)>EPS ) THEN  ! Unique b Try to augment to full H-M
                       spcgr = spcgr(1:1) //'1'// spcgr(2:length) // '1'
                       l_space_group = spcgr_test(spcgr ) ! Test for known space group
                     ELSEIF(ABS(latt(4)-90.0)<EPS .AND. ABS(latt(5)-90.0)<EPS.AND.  &
                        ABS(latt(6)-90.0)>EPS ) THEN  ! Unique c Try to augment to full H-M
                       spcgr = spcgr(1:length) //'11'
                       l_space_group = spcgr_test(spcgr ) ! Test for known space group
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
            IF(l_space_group) THEN
               symm_n = 0
            ELSE                     !, flag error but finish writing
               IF(symm_n>0) THEN 
                  WRITE(iwr, '(a,a,a)') '# spcgr ',spcgr(1:LEN_TRIM(spcgr)), &
                                        ' # Original in CIF is UNKNOWN Used symmetry operations #'
                  spcgr = 'P1'
               ELSE                  !, flag error but finish writing
                  ier_num = -7
                  ier_typ = ER_APPL
               ENDIF
            ENDIF
         ENDIF
         WRITE(iwr, 1100) spcgr(1:LEN_TRIM(spcgr))
      ELSE    ! spcgr is blank
         IF(symm_n == 0) THEN   ! no symmetry operations
            IF(spcgr_no /= 0) THEN
               WRITE(iwr, 1150) spcgr_no
               symm_n = 0
            ELSE
               WRITE(iwr, 1170)
               symm_n = 0
            ENDIF
         ENDIF
      ENDIF
!
      IF(symm_n>0) THEN 
         l_numbers = .FALSE.
         IF(rawline(symm_1)(1:1)=='1' ) THEN
            length = len_trim(rawline(symm_1))
            line = rawline(symm_1)(1:length)
            CALL rem_bl(line, length)
            IF(line(2:2)=='''' .OR. .NOT.(line(2:2)=='/' .OR. line(2:2)=='.')) THEN
               l_numbers = .TRUE.
            ENDIF
         ENDIF
         write_symm: DO i=1, symm_n        ! interpret symmetry lines into matrix form
            IF(l_numbers) THEN
               IF(i<10) THEN 
                 rawline(symm_1+i-1)(1:1) = ' '
               ELSEIF(i<100) THEN 
                 rawline(symm_1+i-1)(1:2) = '  '
               ELSEIF(i<1000) THEN 
                 rawline(symm_1+i-1)(1:3) = '   '
               ENDIF
            ENDIF
            symm_mat(:,:) = 0.0
            iquote1 = INDEX(rawline(symm_1+i-1),'''', .FALSE.)
            IF(iquote1>0) THEN
               iquote2 = INDEX(rawline(symm_1+i-1),'''', .TRUE.)
            ELSE
               iquote2 = LEN_TRIM(rawline(symm_1+i-1))+1
            ENDIF
            length_b= iquote2-iquote1-1
            CALL get_params (rawline(symm_1+i-1) (iquote1+1:iquote2-1), iianz, cspara, lspara, 3, length_b)
            CALL do_cap(cspara(1))
            CALL do_cap(cspara(2))
            CALL do_cap(cspara(3))
            IF(cspara(1)=='X' .AND. cspara(2)=='Y' .AND. cspara(3)=='Z') CYCLE write_symm
            DO j=1,3
               ix = INDEX(cspara(j),'X')
               iy = INDEX(cspara(j),'Y')
               iz = INDEX(cspara(j),'Z')
               IF(ix>0) THEN 
                  symm_mat(j,1) = 1.0
                  cspara(j)(ix:ix) = ' '
                  IF(ix>1.AND.cspara(j)(ix-1:ix-1)=='-') THEN
                     symm_mat(j,1) = -1.0
                     cspara(j)(ix-1:ix-1) = ' '
                  ELSEIF(ix>1.AND.cspara(j)(ix-1:ix-1)=='+') THEN
                     symm_mat(j,1) = +1.0
                     cspara(j)(ix-1:ix-1) = ' '
                  ENDIF
               ENDIF
               IF(iy>0) THEN 
                  symm_mat(j,2) = 1.0
                  cspara(j)(iy:iy) = ' '
                  IF(iy>1.AND.cspara(j)(iy-1:iy-1)=='-') THEN
                     symm_mat(j,2) = -1.0
                     cspara(j)(iy-1:iy-1) = ' '
                  ELSEIF(iy>1.AND.cspara(j)(iy-1:iy-1)=='+') THEN
                     symm_mat(j,2) = +1.0
                     cspara(j)(iy-1:iy-1) = ' '
                  ENDIF
               ENDIF
               IF(iz>0) THEN 
                  symm_mat(j,3) = 1.0
                  cspara(j)(iz:iz) = ' '
                  IF(iz>1.AND.cspara(j)(iz-1:iz-1)=='-') THEN
                     symm_mat(j,3) = -1.0
                     cspara(j)(iz-1:iz-1) = ' '
                  ELSEIF(iz>1.AND.cspara(j)(iz-1:iz-1)=='+') THEN
                     symm_mat(j,3) = +1.0
                     cspara(j)(iz-1:iz-1) = ' '
                  ENDIF
               ENDIF
               IF(cspara(j)/=' ') THEN   ! Additive terms remain
                  ix = INDEX(cspara(j),'.')
                  IF(ix==0) THEN        ! No decimal point
                     lspara(j) = LEN_TRIM(cspara(j))
                     cspara(j)(lspara(j)+1:lspara(j)+1) = '.'
                     lspara(j) = LEN_TRIM(cspara(j))
                  ENDIF
               ELSE
                  cspara(j) = '0.0'
                  lspara(j) = 3
               ENDIF
            ENDDO
            CALL ber_params (3, cspara, lspara, wwerte, 3) 
            symm_mat(1:3,4) = wwerte(1:3)
!
!           Test for unit matrix and omit this, but write all others
!
            symm_mat(1,1) = symm_mat(1,1) - 1.0
            symm_mat(2,2) = symm_mat(2,2) - 1.0
            symm_mat(3,3) = symm_mat(3,3) - 1.0
            IF(MAXVAL(symm_mat)> 0.1 .OR. MINVAL(symm_mat) < -0.1) THEN
               symm_mat(1,1) = symm_mat(1,1) + 1.0
               symm_mat(2,2) = symm_mat(2,2) + 1.0
               symm_mat(3,3) = symm_mat(3,3) + 1.0
               WRITE(iwr, 1180) ((symm_mat(j,k),k=1,4),j=1,3), 1
            ENDIF
         ENDDO write_symm
!
! Currently not implemented as the list of symmetry operations does
! seem to contain centering operations as well
!
! Add centering vectors, just in case, if the list of symmetry operations 
! did not include the centering symmetry operations
! As these are written as generators, they will not blow up the list
! of symmetry operations once the structure is read as cell file
!
!        IF(bravais=='A') THEN
!           WRITE(iwr, 1190) 1.,0.,0.,0.0, 0.,1.,0.,0.5, 0.,0.,1.,0.5, 1
!        ELSEIF(bravais=='B') THEN
!           WRITE(iwr, 1190) 1.,0.,0.,0.5, 0.,1.,0.,0.0, 0.,0.,1.,0.5, 1
!        ELSEIF(bravais=='C') THEN
!           WRITE(iwr, 1190) 1.,0.,0.,0.0, 0.,1.,0.,0.5, 0.,0.,1.,0.5, 1
!        ELSEIF(bravais=='I') THEN
!           WRITE(iwr, 1190) 1.,0.,0.,0.5, 0.,1.,0.,0.5, 0.,0.,1.,0.5, 1
!        ELSEIF(bravais=='F') THEN
!           WRITE(iwr, 1190) 1.,0.,0.,0.0, 0.,1.,0.,0.5, 0.,0.,1.,0.5, 1
!           WRITE(iwr, 1190) 1.,0.,0.,0.5, 0.,1.,0.,0.0, 0.,0.,1.,0.5, 1
!           WRITE(iwr, 1190) 1.,0.,0.,0.5, 0.,1.,0.,0.5, 0.,0.,1.,0.0, 1
!        ELSEIF(bravais=='R') THEN
!           WRITE(iwr, 1190) 1.,0.,0.,2./3., 0.,1.,0.,1./3., 0.,0.,1.,1./3., 1
!           WRITE(iwr, 1190) 1.,0.,0.,1./3., 0.,1.,0.,2./3., 0.,0.,1.,2./3., 1
!        ELSEIF(bravais=='O') THEN
!           WRITE(iwr, 1190) 1.,0.,0.,2./3., 0.,1.,0.,1./3., 0.,0.,1.,2./3., 1
!           WRITE(iwr, 1190) 1.,0.,0.,1./3., 0.,1.,0.,2./3., 0.,0.,1.,1./3., 1
!        ENDIF
      ENDIF
      WRITE(iwr, 1200) latt
      WRITE(iwr, 1300)
      TAIL => HEAD
      TEMP => HEAD
      DO WHILE (ASSOCIATED(TAIL))
         WRITE(iwr,1400) TAIL%at_name,TAIL%at_pos,TAIL%at_bvalue, TAIL%at_occ 
         TAIL => TAIL%next
         DEALLOCATE(TEMP)       ! Clean up the memory structure
         TEMP => TAIL
      ENDDO
      NULLIFY(HEAD)
      NULLIFY(TEMP)
      NULLIFY(TAIL)
1000 FORMAT('title ',a)
1100 FORMAT('spcgr ',a)
1150 FORMAT('spcgr ',i5)
1170 FORMAT('spcgr  P1')
1190 FORMAT('gene  ',3(3(f5.1,', '),f12.9,','), I3)
1180 FORMAT('symm  ',3(3(f5.1,', '),f12.9,','), I3)
1200 FORMAT('cell  ',5(f12.5,', '),f12.5)
1300 FORMAT('atoms x,',12x,'y,',12x,'z,',12x,'Biso,', 4x,'Property,', &
                 2x,'MoleNo,  MoleAt, Occ')
1400 FORMAT(a4, 4(F12.8,', '),'1,      0,       0,      ',F8.6  )
!
      CLOSE(iwr)
!
         data_i = data_i + 1      ! We wrote a data section, increment counter
         title = newtitle         ! title will be written to the next file
      ENDDO data_entries
!
! clean up arrays
!
      DEALLOCATE(rawline)
!
      END SUBROUTINE cif2discus
!
!*******************************************************************************
!
subroutine lammps2discus(ianz, cpara, lpara, MAXW, &
                            NOPTIONAL, opara, lopara, lpresent, O_ATOM, ofile) 
!-                                                                      
!     converts a LAMMPS file to DISCUS                   
!  Presently a very basic routine for a triclinic dump.
!+                                                                      
!                                                                       
!use tensors_mod
use build_name_mod
!use wink_mod
!use ber_params_mod
!use blanks_mod
use get_params_mod
!use lib_length
use precision_mod
!use string_convert_mod
use support_mod
use trig_degree_mod
!
implicit none 
!                                                                       
integer                                , intent(inout) :: ianz 
integer                                , intent(in)    :: MAXW 
character (LEN=*), DIMENSION(1:MAXW)   , intent(inout) :: cpara
integer          , DIMENSION(1:MAXW)   , intent(inout) :: lpara
integer                                , intent(in)    :: NOPTIONAL
character(LEN=*) , dimension(NOPTIONAL), intent(in)    :: opara   !Optional parameter strings returned
integer          , dimension(NOPTIONAL), intent(in)    :: lopara  !Lenght opt. para name returned
logical          , dimension(NOPTIONAL), intent(in)    :: lpresent!opt. para is present
integer                                , intent(in)    :: O_ATOM
character(LEN=*)                       , intent(out)   :: ofile 
!                                                                       
integer, parameter :: ird = 34            ! I/O channel read
integer, parameter :: iwr = 35            ! I/O channel write
integer, parameter :: MAXP = 100          ! Max number of parameters on atoms:[]
!
character(len=PREC_STRING) :: infile      ! Input file 
character(len=PREC_STRING) :: line        ! Dummy line 
character(len=PREC_STRING) :: string      ! Dummy line 
character(len=PREC_STRING) :: title       ! Crystal structure title line
integer :: i, j                           ! Dummy counters
integer :: ios                            ! I/O status flag
integer :: natoms                         ! Grand number of atoms
integer :: nr                             ! Atom number
integer :: id                             ! Atom type
integer :: length                         ! Dummy length
real(kind=PREC_DP), dimension(3)    :: xyz  ! Atom coordinates bounding box
real(kind=PREC_DP), DIMENSION(MAXW) :: werte
real(kind=PREC_DP) :: xlo_bound, xhi_bound, xy  ! Bounding box parameters 
real(kind=PREC_DP) :: ylo_bound, yhi_bound, xz
real(kind=PREC_DP) :: zlo_bound, zhi_bound, yz
real(kind=PREC_DP) :: lx, ly, lz                ! Bounding box lengths
real(kind=PREC_DP), DIMENSION(3,2 ) :: xyz_lh   ! LAMMPS xlow xhigh etc
real(kind=PREC_DP), DIMENSION(6)    :: lat      ! lattice parametersetc
!
character(LEN=4) , dimension(:), allocatable :: catom   !Optional parameter strings returned
integer          , dimension(:), allocatable :: latom   !Lenght opt. para name returned
!
!                                                                       
!     Create input / output file name
!
werte = 0.0D0
CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
IF (ier_num.ne.0) THEN 
   RETURN 
ENDIF 
infile = cpara (1) 
i = index (infile, '.',.true.)                  ! find last occurence of '.'
IF (i.eq.0) THEN 
   infile = cpara(1)(1:lpara(1)) //'.out' 
   ofile  = cpara(1)(1:lpara(1)) //'.stru' 
   title  = cpara(1)(1:lpara(1))
ELSE 
   IF(    cpara(1)(lpara(1)-3:lpara(1)) == '.out') THEN
      ofile  = cpara (1) (1:lpara(1)-3) //'stru' 
   ELSEIF(cpara(1)(lpara(1)-3:lpara(1)) == '.OUT') THEN
      ofile  = cpara (1) (1:lpara(1)-3) //'stru' 
   ELSE
      ofile  = cpara (1) (1:i) //'stru'
   ENDIF
ENDIF 
!
CALL oeffne(ird, infile, 'old') 
IF (ier_num.ne.0) THEN 
   RETURN 
ENDIF 
!
!  Read header section 
!
header: do
   read(ird,'(a)', iostat=ios) line
   if(is_iostat_end(ios)) then
      ier_num = 6
      ier_typ = ER_IO
      ier_msg(1) = 'LAMMPS header ended unexpectedly'
      close(ird)
      return
   endif
   if(line(1:5)=='ITEM:') then
      if(line(7:14)=='TIMESTEP') then
         read(ird,'(a)', iostat=ios) line
         title = title(1:len_trim(title)) // ' LAMMPS Timestep ' // line(1:len_trim(line))
      elseif(line(7:21)=='NUMBER OF ATOMS') then
         read(ird,'(a)', iostat=ios) line
         read(line,*) natoms
      elseif(line(7:16)=='BOX BOUNDS') then
         if(line(18:25)=='xy xz yz') then
            read(ird,'(a)', iostat=ios) line
            read(line,*) xlo_bound, xhi_bound, xy
            read(ird,'(a)', iostat=ios) line
            read(line,*) ylo_bound, yhi_bound, xz
            read(ird,'(a)', iostat=ios) line
            read(line,*) zlo_bound, zhi_bound, yz
!
            xyz_lh(1,1) = xlo_bound - min(0.0D0,xy,xz,xy+xz)
            xyz_lh(1,2) = xhi_bound - max(0.0D0,xy,xz,xy+xz)
            xyz_lh(2,1) = ylo_bound - min(0.0D0,yz)
            xyz_lh(2,2) = yhi_bound - max(0.0D0,yz)
            xyz_lh(3,1) = zlo_bound
            xyz_lh(3,2) = zhi_bound
            lx = xyz_lh(1,2) - xyz_lh(1,1)
            ly = xyz_lh(2,2) - xyz_lh(2,1)
            lz = xyz_lh(3,2) - xyz_lh(3,1)
            lat(1) = lx
            lat(2) = sqrt( ly**2 + xy**2 )
            lat(3) = sqrt( lz**2 + xz**2 + yz**2 )
            lat(4) = acosd( (xy*xz + lx*yz)/lat(2)/lat(3))
            lat(5) = acosd( xz/lat(3) )
            lat(6) = acosd( xy/lat(2) )
         endif
      elseif(line(7:11)=='ATOMS') then
         exit header                     ! Last header line, start to read atoms
      endif
   endif
enddo header
!
!  set up atom names either dummy Q1... or from optional parameters
!
allocate(catom(1:MAXP))
allocate(latom(1:MAXP))
if(lpresent(O_ATOM)) then
   length = lopara(O_ATOM)
   string = opara(O_atom)(2:lopara(O_ATOM)-1)
   call get_params(string, ianz, catom, latom, MAXP, length)
else
   do i=1,MAXP
      write(catom(i),'(a,i3.3)') 'Q',i
   enddo
endif
!
open(IWR, file=ofile, status='unknown')
write(IWR, '(a,a)') 'title ', title(1:len_trim(title))
write(IWR, '(a)  ') 'spcgr P1'
write(IWR, '(a, 5(f10.6,a),f10.6)') 'cell ',(lat(i),',',i=1,5), lat(6)
write(IWR, '(a,a)') 'scat '
write(IWR, '(a,a)') 'adp  '
write(IWR, '(a,a)') 'occ  '
write(IWR, '(a)') 'atoms'
!
!  Loop over atoms
!
atoms: do j=1, natoms
   read(ird,'(a)', iostat=ios) line
   if(is_iostat_end(ios)) then
      ier_num = 6
      ier_typ = ER_IO
      ier_msg(1) = 'LAMMPS atoms ended unexpectedly'
      close(IRD)
      close(IWR)
      deallocate(catom)
      deallocate(latom)
      return
   endif
   read(line,*) nr, id, xyz
   write(IWR, '(a4, 2x,3(f9.6,a),f6.2)') catom(id), (xyz(i),',',i=1,3), 1.0
enddo atoms
close(IWR)
close(IRD)
deallocate(catom)
deallocate(latom)
!
end subroutine lammps2discus
!
!*******************************************************************************
!
      SUBROUTINE test_atom_name(lsymbol, at_name)
!
      USE element_data_mod
      USE charact_mod
      USE blanks_mod
      USE errlist_mod
USE lib_errlist_func
      USE string_convert_mod
      IMPLICIT NONE
!
      LOGICAL,          INTENT(IN)    :: lsymbol
      CHARACTER(LEN=*), INTENT(INOUT) :: at_name
!
      CHARACTER (LEN=2), DIMENSION(-8:8) :: ions
      CHARACTER (LEN=2), DIMENSION(-8:8) :: names
      INTEGER :: length, i, j, charge, el_number
!
      DATA ions / '-8', '-7', '-6', '-5', '-4', '-3', '-2', '-1', '  ',  &
                  '+1', '+2', '+3', '+4', '+5', '+6', '+7', '+8' /
      DATA names/ '8-', '7-', '6-', '5-', '4-', '3-', '2-', '1-', '  ',  &
                  '1+', '2+', '3+', '4+', '5+', '6+', '7+', '8+' /
!
      length = LEN_TRIM(at_name)   
      charge = 0
      CALL do_cap(at_name)
!
      IF(lsymbol) THEN             ! Atom symbols, check charge
         charge = MAX(INDEX(at_name,'-'), INDEX(at_name,'+'))
         IF(charge>0) THEN         ! A charge is given, check order and entry
            IF(charge==length-1) THEN
               DO j=-8,8
                  IF(at_name(charge:charge+1)==ions(j)) at_name(charge:charge+1) = names(j)
               ENDDO
            ELSEIF(charge==length) THEN
               DO j=-8,8
                  IF(at_name(charge-1:charge)==names(j)) at_name(charge-1:charge) = names(j)
               ENDDO
            ELSEIF(IACHAR(at_name(charge-1:charge-1))<=zero .OR.  nine>=IACHAR(at_name(charge-1:charge-1))) THEN
               at_name(charge:length) = ' '
            ENDIF
            CALL rem_bl(at_name,length)
            CALL symbf(at_name, el_number)
            IF(el_number > 0) RETURN
!                                  ! Element not found, try without charge
            CALL symbf(at_name(1:charge-1), el_number)
            IF(el_number > 0) THEN
               ier_num = +5
               ier_typ = ER_APPL
               ier_msg(1) = 'Atom '//at_name//' renamed to '//at_name(1:charge-1)
               CALL errlist
               at_name = at_name(1:charge-1)
               RETURN
            ENDIF
            ier_num = +6
            ier_typ = ER_APPL
            ier_msg(1) = 'Atom '//at_name//' unknown '
            CALL errlist
         ENDIF
      ELSE                         ! Atom labels, ignore anything but characters
         DO i=1,length
            IF(IACHAR(at_name(i:i))<AA .OR. IACHAR(at_name(i:i))>zz) at_name(i:i) = ' '
         ENDDO
         CALL rem_bl(at_name,length)
      ENDIF
!
!
      END SUBROUTINE test_atom_name
!
!*******************************************************************************
!
SUBROUTINE test_file ( strucfile, natoms, ntypes, n_mole, n_type, &
                       n_atom, n_cells, init, l_cell)
!
!     Determines the number of atoms and atom types in strucfile
!
use atom_line_mod
!
use blanks_mod
USE ber_params_mod
USE charact_mod
USE get_params_mod
USE lib_length
USE precision_mod
USE str_comp_mod
USE string_convert_mod
USE support_mod
!
IMPLICIT NONE
!
CHARACTER (LEN=*)    , INTENT(IN)    :: strucfile
INTEGER              , INTENT(INOUT) :: natoms
INTEGER              , INTENT(INOUT) :: ntypes
INTEGER              , INTENT(INOUT) :: n_mole 
INTEGER              , INTENT(INOUT) :: n_type 
INTEGER              , INTENT(INOUT) :: n_atom 
integer, dimension(3), intent(out)   :: n_cells
INTEGER              , INTENT(IN)    :: init
LOGICAL              , INTENT(IN)    :: l_cell         ! If true this is a 'cell' > all atoms differ in type
!
INTEGER, PARAMETER                    :: MAXW = 16 
CHARACTER(LEN=PREC_STRING), DIMENSION(MAXW)  :: cpara (MAXW) 
INTEGER            , DIMENSION(MAXW)  :: lpara (MAXW) 
REAL(KIND=PREC_DP) , DIMENSION(MAXW)  :: werte (MAXW) 
!
REAL, PARAMETER                       :: eps = 1e-6
CHARACTER (LEN=PREC_STRING)                  :: line
CHARACTER (LEN=PREC_STRING)                  :: zeile
CHARACTER (LEN=  20)                  :: bef
CHARACTER (LEN=   4), DIMENSION(PREC_STRING), SAVE :: names
REAL                , DIMENSION(PREC_STRING), SAVE :: bvals
REAL                , DIMENSION(PREC_STRING), SAVE :: occs
INTEGER                               :: ios
INTEGER                               :: i
INTEGER                               :: iblk
INTEGER                               :: ianz   ! no of arguments
INTEGER                               :: laenge ! length of input line
INTEGER                               :: lp     ! length of parameter string
INTEGER                               :: nscattypes ! no of SCAT arguments 
INTEGER                               :: nadptypes  ! no of ADP  arguments 
INTEGER                               :: nocctypes  ! no of ADP  arguments 
INTEGER                               :: indxt      ! Pos of a TAB in input
INTEGER                               :: indxb      ! Pos of a BLANK in input
INTEGER                               :: lbef       ! Length of command string
INTEGER                               :: iflag      ! FLAG FOR NMDE 
INTEGER                               :: imole      ! Mole number on atom line
INTEGER                               :: inatom     ! Atom in Mole number on atom line
LOGICAL                               :: in_mole    ! Currently within a molecule
LOGICAL                               :: l_type     ! RFound molecule type command
integer :: nlines ! number of line sin input file
LOGICAL                               :: new
REAL                                  :: xc,yc,zc,bval
real(kind=PREC_SP), dimension(3,2)   :: ccdim     ! Crystal dimensions 
REAL                                 :: occ
integer, dimension(9) :: ncell_val
!
LOGICAL           :: IS_IOSTAT_END
!
ncell_val  = 0
natoms     = 0
nscattypes = 0
nadptypes  = 0
nocctypes  = 0
!n_mole     = 0
ntypes     = 0
iblk = 0
ccdim = 0.0
IF ( init == -1 ) THEN
   names(:)  = ' '
   bvals(:)  = 0.0
   occs(:)   = 1.0
   ntypes    = 0
   n_mole     = 0
   n_type     = 0
   n_atom     = 0
ENDIF
in_mole = .false.
!
call atom_get_size(strucfile, nlines)
call atom_alloc(nlines)
!
CALL oeffne ( 99, strucfile, 'old')
IF ( ier_num /= 0) THEN
   ier_msg(3) = strucfile(MAX(1,LEN_TRIM(strucfile)-LEN(ier_msg)):LEN_TRIM(strucfile))
   CLOSE ( 99 )
   RETURN
ENDIF
!
header: DO
   READ (99,'(a)', IOSTAT=ios) line
   IF(ios /= 0) THEN
      ier_num = -6
      ier_typ = ER_IO
      CLOSE ( 99 )
      ier_msg(3) = strucfile(MAX(1,LEN_TRIM(strucfile)-LEN(ier_msg)):LEN_TRIM(strucfile))
      call atom_dealloc
      RETURN
   ENDIF
   IF (line == ' '.OR.line (1:1)  == '#'.OR. line(1:1) == '!' .OR. &
       line == CHAR (13) )  CYCLE header
   laenge = LEN_TRIM(LINE)
   call tab2blank(line, laenge)
   iblk = INDEX(line, ' ')
   CALL do_cap (line(1:iblk))
!        laenge = len_str(line)
   IF ( laenge .gt. 4 ) THEN
      zeile = line(iblk+1:laenge)
      lp    = laenge - 4
   ELSE
      zeile = ' '
      lp    = 1
   ENDIF
   IF (line(1:4) == 'SCAT' ) THEN 
       CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
       IF (ier_num.eq.0) THEN 
          DO i = 1,ianz
              names(nscattypes+i) = cpara(i)(1:lpara(i))
          ENDDO
          nscattypes = nscattypes + ianz
       ELSE
          ier_num = -111
          ier_typ = ER_APPL
          ier_msg(3) = strucfile(MAX(1,LEN_TRIM(strucfile)-LEN(ier_msg)):LEN_TRIM(strucfile))
          CLOSE(99)
      call atom_dealloc
          RETURN
       ENDIF
   ELSEIF (line(1:3) == 'ADP' ) THEN 
       CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
       IF (ier_num.eq.0) THEN 
          CALL ber_params (ianz, cpara, lpara, werte, maxw) 
          IF (ier_num.eq.0) THEN 
             DO i = 1,ianz
                 bvals(nadptypes+i) = werte(i)
             ENDDO
             nadptypes = nadptypes + ianz
          ELSE
             ier_num = -112
             ier_typ = ER_APPL
             ier_msg(3) = strucfile(MAX(1,LEN_TRIM(strucfile)-LEN(ier_msg)):LEN_TRIM(strucfile))
             CLOSE(99)
      call atom_dealloc
             RETURN
          ENDIF
       ELSE
          ier_num = -149
          ier_typ = ER_APPL
          ier_msg(3) = strucfile(MAX(1,LEN_TRIM(strucfile)-LEN(ier_msg)):LEN_TRIM(strucfile))
          CLOSE(99)
      call atom_dealloc
          RETURN
       ENDIF
   ELSEIF (line(1:3) == 'OCC' ) THEN 
       CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
       IF (ier_num.eq.0) THEN 
          CALL ber_params (ianz, cpara, lpara, werte, maxw) 
          IF (ier_num.eq.0) THEN 
             DO i = 1,ianz
                 IF(werte(i)<0.0 .OR. 1.0<werte(i)) THEN
                    ier_num = -150
                    ier_typ = ER_APPL
                    ier_msg(3) = strucfile(MAX(1,LEN_TRIM(strucfile)-LEN(ier_msg)):LEN_TRIM(strucfile))
                    CLOSE(99)
                    RETURN
                 ENDIF
                 occs (nocctypes+i) = werte(i)
             ENDDO
             nocctypes = nocctypes + ianz
          ELSE
             ier_num = -149
             ier_typ = ER_APPL
             ier_msg(3) = strucfile(MAX(1,LEN_TRIM(strucfile)-LEN(ier_msg)):LEN_TRIM(strucfile))
             CLOSE(99)
      call atom_dealloc
             RETURN
          ENDIF
       ELSE
          ier_num = -149
          ier_typ = ER_APPL
          ier_msg(3) = strucfile(MAX(1,LEN_TRIM(strucfile)-LEN(ier_msg)):LEN_TRIM(strucfile))
          CLOSE(99)
      call atom_dealloc
          RETURN
       ENDIF
   ELSEIF (line(1:4) == 'CELL' ) THEN 
       CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
       IF (ier_num.eq.0 .AND. ianz == 6) THEN 
          CALL ber_params (ianz, cpara, lpara, werte, maxw) 
          IF (ier_num /= 0) THEN 
             ier_num = -48
             ier_typ = ER_APPL
             ier_msg(3) = strucfile(MAX(1,LEN_TRIM(strucfile)-LEN(ier_msg)):LEN_TRIM(strucfile))
             CLOSE(99)
      call atom_dealloc
             RETURN
          ENDIF
       ELSE
          READ(zeile,*,IOSTAT=ios) (werte(i),i=1,6)
          IF(ios /=0 .OR. is_nan(REAL(werte(1),KIND=PREC_SP)) .OR. is_nan(REAL(werte(2),KIND=PREC_SP)) &
                     .OR. is_nan(REAL(werte(3),KIND=PREC_SP)) .OR. is_nan(REAL(werte(4),KIND=PREC_SP)) &
                     .OR. is_nan(REAL(werte(5),KIND=PREC_SP)) .OR. is_nan(REAL(werte(6),KIND=PREC_SP))) THEN
             ier_num = -48
             ier_typ = ER_APPL
             ier_msg(3) = strucfile(MAX(1,LEN_TRIM(strucfile)-LEN(ier_msg)):LEN_TRIM(strucfile))
             CLOSE(99)
      call atom_dealloc
             RETURN
          ENDIF
       ENDIF
   ELSEIF (line(1:5) == 'NCELL' ) THEN 
       CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
       IF (ier_num.eq.0) THEN 
          werte = 0.0
          CALL ber_params (ianz, cpara, lpara, werte, maxw) 
          IF (ier_num.eq.0) THEN 
             DO i = 1,ianz
                 IF(werte(i)<0.0 ) THEN
                    ier_num = -164
                    ier_typ = ER_APPL
                    ier_msg(3) = strucfile(MAX(1,LEN_TRIM(strucfile)-LEN(ier_msg)):LEN_TRIM(strucfile))
                    CLOSE(99)
                    call atom_dealloc
                    RETURN
                 ENDIF
             ENDDO
             IF(ianz>=5) THEN
                natoms = werte(5)
             ENDIF
             n_cells = werte(1:3)
             ntypes = werte(6)
             n_mole = werte(7)
             n_type = werte(8)
             n_atom = werte(9)
             ncell_val = nint(werte(1:9))
          ELSE
             ier_num = -163
             ier_typ = ER_APPL
             ier_msg(3) = strucfile(MAX(1,LEN_TRIM(strucfile)-LEN(ier_msg)):LEN_TRIM(strucfile))
             CLOSE(99)
             call atom_dealloc
             RETURN
          ENDIF
       ELSE
          ier_num = -163
          ier_typ = ER_APPL
          ier_msg(3) = strucfile(MAX(1,LEN_TRIM(strucfile)-LEN(ier_msg)):LEN_TRIM(strucfile))
          CLOSE(99)
          call atom_dealloc
          RETURN
       ENDIF
   ENDIF
   IF(line(1:4) == 'ATOM') THEN 
      call atom_line_inter(zeile, lp)
      EXIT header
   ENDIF
ENDDO header
!
IF (nscattypes /= nadptypes ) THEN
   ier_num = -115
   ier_typ = ER_APPL
   ier_msg(3) = strucfile(MAX(1,LEN_TRIM(strucfile)-LEN(ier_msg)):LEN_TRIM(strucfile))
   CLOSE(99)
   call atom_dealloc
   RETURN
ENDIF
!
ntypes = MAX(ntypes,nscattypes)
!
IF(natoms>0) THEN
   IF(ntypes>0) THEN
      IF(n_mole>0) THEN
         n_type = MAX(n_type,1)  ! Ensure values are NOT zero
         n_atom = MAX(n_atom,1)
      ENDIF
!     CLOSE(99)
!     RETURN
   ENDIF
ENDIF
if(ncell_val(5) > 0) then
   natoms = 0    ! temporarily set back to zero, as vals are in ncell_val
   ntypes = 0
   n_mole = 0
   n_type = 0
   n_atom = 0
endif
!
l_type = .FALSE.
main: DO
   READ (99,'(a)', IOSTAT=ios) line
   IF( IS_IOSTAT_END(ios) ) EXIT main
   IF(line == ' '.OR.line (1:1)  == '#'.OR. line(1:1) == '!' .OR. &
      line == CHAR (13) )  CYCLE main
   laenge = len_str(line)
   call tab2blank(line, laenge)
!
   bef   = '    '
   indxt = INDEX (line, tab)       ! find a tabulator
   IF(indxt==0) indxt = laenge + 1
   indxb = index (line, ' ')       ! find a blank
   IF(indxb==0) indxb = laenge + 1
   indxb = MIN(indxb,indxt)
!
   lbef = min (indxb - 1, 8)
   bef  = line (1:lbef)
   CALL do_cap (line(1:laenge))
   IF(line(1:1)=='#' .OR. line(1:1)=='!') CYCLE main
!
   ismole: IF(str_comp(line, 'MOLECULE', 3, lbef, 8) .or. &
              str_comp(line, 'DOMAIN'  , 3, lbef, 6) .or. &
              str_comp(line, 'OBJECT'  , 3, lbef, 6)     ) THEN
      IF ( indxb+1 >= laenge) THEN   ! No parameter => start
         IF ( .not. in_mole) THEN
            in_mole = .true.
            n_mole  = n_mole + 1
            l_type  = .false.
         ENDIF
      ELSEIF ( str_comp(line(indxb+1: laenge), 'END',3, laenge-indxb,3)) THEN
         IF ( in_mole) THEN
            in_mole = .false.
            IF(.not.l_type) THEN
               n_type = n_type + 1
            ENDIF
            l_type  = .false.
         ENDIF
      ELSEIF ( str_comp(line(indxb+1: laenge), 'TYPE',3, laenge-indxb,4)) THEN
         zeile  = line(indxb+1: laenge)
         laenge = laenge-indxb
         CALL get_params (zeile, ianz, cpara, lpara, maxw, laenge)
         cpara(1) = '0' 
         lpara(1) = 1
         CALL ber_params (ianz, cpara, lpara, werte, maxw)
         n_type = MAX(n_type, NINT(werte(2)))
         l_type = .true.
!          ELSE
      ENDIF
   ELSE ismole
      if(at_style==-1) then
         call atom_line_get_style(line, lbef+1, laenge, MAXW, werte)
      endif
      iflag  = 0
      imole  = 0
      inatom = 0
      ios = 0
      CALL read_atom_line (line, lbef+1, laenge, natoms, MAXW, werte)! , &
!                           AT_MAXP, at_ianz, at_param, at_init)                                          
      xc     = werte(1)
      yc     = werte(2)
      zc     = werte(3)
      bval   = werte(4)
      iflag  = NINT(werte(5))
      imole  = NINT(werte(6))
      inatom = NINT(werte(7))
      occ    =      werte(8)
      IF(occ     <0.0 .OR. 1.0<occ     ) THEN
         ier_num = -150
         ier_typ = ER_APPL
         ier_msg(1) = line(1:46)
         WRITE(ier_msg(2),'(a,i8)') 'Atom nr. ', natoms + 1
         ier_msg(3) = strucfile(MAX(1,LEN_TRIM(strucfile)-LEN(ier_msg)):LEN_TRIM(strucfile))
         CLOSE(99)
         call atom_dealloc
         RETURN
      ENDIF
         n_mole = MAX(n_mole, imole)
         n_atom = MAX(n_atom, inatom)
      IF(is_nan(xc) .OR. is_nan(yc) .OR. is_nan(zc) .OR. is_nan(bval)) THEN
         ios = -1
      ENDIF
      IF(     is_nan(REAL(werte(10),KIND=PREC_SP)) &
         .OR. is_nan(REAL(werte(11),KIND=PREC_SP)) &
         .OR. is_nan(REAL(werte(12),KIND=PREC_SP)) ) THEN
         ios = -1
      ENDIF
!
      isatom:    IF ( ios == 0 ) THEN
         ccdim(1,1) = min(ccdim(1,1), xc)
         ccdim(1,2) = max(ccdim(1,2), xc)
         ccdim(2,1) = min(ccdim(2,1), yc)
         ccdim(2,2) = max(ccdim(2,2), yc)
         ccdim(3,1) = min(ccdim(3,1), zc)
         ccdim(3,2) = max(ccdim(3,2), zc)
         natoms = natoms + 1
         IF ( in_mole ) THEN
            n_atom = n_atom + 1
         ENDIF
         new = .true.
         if(.not.l_cell) then
         types: DO i=1,ntypes
               if(line(1:lbef) == names(i) .and.                 &
                  abs(abs(bval)-abs(bvals(i))) < eps  .and.      &
                  abs(abs(occ )-abs( occs(i))) < eps        ) then
                  new = .false.
                  exit types
               endif
!           IF ( line(1:lbef) == names(i) ) THEN
!              IF ( l_cell ) THEN
!                 new = .false.
!                 EXIT types
!              ELSEIF(abs(abs(bval)-abs(bvals(i))) < eps .AND.   &
!                     abs(abs(occ )-abs( occs(i))) < eps ) THEN
!                 new = .false.
!                 EXIT types
!              ENDIF
!           ENDIF
         ENDDO types
         endif
         IF ( new ) THEN
            ntypes = ntypes + 1
            names(ntypes) = line(1:lbef)
            bvals(ntypes) = bval
            occs(ntypes) = occ 
         ENDIF
      ELSE isatom
         ier_num = -49
         ier_typ = ER_APPL
         ier_msg(1) = line(1:46)
         WRITE(ier_msg(2),'(a,i8)') 'Atom nr. ', natoms + 1
         ier_msg(3) = strucfile(MAX(1,LEN_TRIM(strucfile)-LEN(ier_msg)):LEN_TRIM(strucfile))
         CLOSE(99)
         call atom_dealloc
         RETURN
      ENDIF isatom
   ENDIF ismole
ENDDO main
!
CLOSE (99)
!
IF(n_mole>0) THEN
   n_type = MAX(n_type,1)  ! Ensure values are NOT zero
   n_atom = MAX(n_atom,1)
ENDIF
if(ncell_val(5)>0) then
   natoms = max(natoms, ncell_val(5))
   ntypes = max(ntypes, ncell_val(6))
   n_mole = max(n_mole, ncell_val(7))
   n_type = max(n_type, ncell_val(8))
   n_atom = max(n_atom, ncell_val(9))
endif
!
n_cells     = int(ccdim(:,2)-ccdim(:,1)) + 1
!
END SUBROUTINE test_file
!
!*******************************************************************************
!
SUBROUTINE test_identical(l_identical, r_identical)
!
USE crystal_mod
USE metric_mod
USE errlist_mod
use precision_mod
IMPLICIT NONE
!
LOGICAL         ,                  INTENT(IN) :: l_identical
REAL            ,                  INTENT(IN) :: r_identical
!
LOGICAL, PARAMETER :: LSPACE = .TRUE.
REAL    :: eps
REAL(kind=PREC_DP), DIMENSION(3) :: u,v
INTEGER :: i, j
!
eps = 1.0E-5
IF(l_identical) eps = r_identical
!
main: DO i=1, cr_natoms-1
   u= cr_pos(:,i)
   DO j=i+1,cr_natoms
      v= cr_pos(:,j)
      IF(do_blen(lspace, u,v)< EPS) THEN
!     IF(ABS(cr_pos(1,i)-cr_pos(1,j)) < eps .AND.  &
!        ABS(cr_pos(2,i)-cr_pos(2,j)) < eps .AND.  &
!        ABS(cr_pos(3,i)-cr_pos(3,j)) < eps ) THEN
         ier_num = -141
         ier_typ = ER_APPL
         WRITE(ier_msg(1),'(a,i6, a, i6,a)') 'Atoms ',i,',',j, &
                          ' are at identical positions' 
         ier_msg(2) = 'Atoms might be separated by integer unit cells'
         IF(l_identical) THEN
            ier_msg(3) = 'If intended, decrease: radius:value '
         ELSE
            ier_msg(3) = 'If intended, use: identical:tolerate'
         EXIT main
         ENDIF
      ENDIF
   ENDDO
ENDDO main
END SUBROUTINE test_identical
!
!*******************************************************************************
!
SUBROUTINE set_spcgr(line,length)
!
USE crystal_mod
USE spcgr_apply
USE gen_add_mod
USE sym_add_mod
USE wyckoff_mod
USE errlist_mod
USE get_params_mod
USE precision_mod
USE take_param_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: length
!
INTEGER, PARAMETER                   ::MAXW = 2
!
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(MAXW) :: cpara
INTEGER            , DIMENSION(MAXW) :: lpara
REAL(KIND=PREC_DP) , DIMENSION(MAXW) :: werte
INTEGER :: ianz
!
INTEGER, PARAMETER :: NOPTIONAL = 1
INTEGER, PARAMETER :: O_SETTING = 1
CHARACTER(LEN=   7)       , DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 0 ! Number of values to calculate 
!                                                                       
!
DATA oname  / 'setting' /
DATA loname /  7        /
!
opara  =  (/ 'abc' /)   ! Always provide fresh default values
lopara =  (/  3   /)
owerte =  (/  0.0 /)
!
!
CALL get_params (line, ianz, cpara, lpara, maxw, length)
IF(ianz>=1) THEN                   ! At least one parameter
   cr_spcgr= cpara(1)(1:lpara(1))  ! Set space group name
   IF(ianz==2) THEN
      spcgr_para = NINT(werte(2))  ! Set origin choice parameter
      spcgr_ianz = 2
   ELSE
      spcgr_ianz = 1
   ENDIF
   werte (1) = spcgr_para 
!  Optionally get setting
   CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                     oname, loname, opara, lopara, lpresent, owerte)
   cr_set = opara(O_SETTING)(1:MIN(3,lopara(O_SETTING)))
   CALL spcgr_no (spcgr_ianz, MAXW, werte) 
   IF(ier_num/=0) RETURN
   spc_n     = 1     ! reset spcgr
   gen_add_n = 0     ! No additional generators
   sym_add_n = 0     ! No additional symmetry operations
   CALL get_symmetry_matrices
   IF(cr_syst/=4 .AND. cr_iset/=1) THEN !non-orthorhombic and nonstandard setting
      ier_num = -160
      ier_typ = ER_APPL
      ier_msg(1) = 'Non-standard settings are implemented only'
      ier_msg(2) = 'for orthorhombic space groups '
      RETURN
   ENDIF
ELSE
   ier_num = -6
   ier_typ = ER_COMM
ENDIF
!
END SUBROUTINE set_spcgr
!
!*******************************************************************************
!
SUBROUTINE test_mole_gap
!
USE molecule_mod
!
USE errlist_mod
!
INTEGER :: i,j
!
main: DO i=1,MIN(mole_num_atom,UBOUND(mole_cont,1))
   IF(mole_cont(i) == 0) THEN
      ier_msg(1) = 'Atoms are missing in a molecule'
      second: DO j=1,mole_num_mole
         IF(i<mole_off(j)+mole_len(j)) THEN
            WRITE(ier_msg(1),'(A,i7)') 'Atoms are missing in molecule', j
            EXIT second
         ENDIF
      ENDDO second
      ier_num = -156
      ier_typ = ER_APPL
      ier_msg(2) = 'Check the MoleAt column in the structure file'
      ier_msg(3) = 'Was structure saved with atoms deselected?'
      EXIT main
   ENDIF
ENDDO main
!
IF(ier_num/=0) CALL rese_cr
!
END SUBROUTINE test_mole_gap
!
!*******************************************************************************
!
SUBROUTINE readcell_mole(strucfile, l_identical, r_identical)
!
USE crystal_mod
USE discus_save_mod
USE molecule_mod
USE prop_para_func
USE read_internal_mod
USE spcgr_apply
USE save_menu, ONLY: save_internal, save_store_setting, save_restore_setting, save_default_setting, save_struc, save_show
use wyckoff_mod
USE lib_errlist_func
USE precision_mod
!
IMPLICIT none 
!
!                                                                       
CHARACTER ( LEN=* ), INTENT(OUT) :: strucfile 
LOGICAL            , INTENT(IN)  :: l_identical
REAL               , INTENT(IN)  :: r_identical
!
CHARACTER(LEN=200) :: tempfile
CHARACTER(LEN=PREC_STRING) :: line
INTEGER :: length
INTEGER :: im, j, iat
LOGICAL :: lout
INTEGER, DIMENSION(3) :: n_unit_cells  ! local copy to survive readstru 
REAL, DIMENSION(3) :: vec     ! position of first atom in a molecule
REAL, DIMENSION(3) :: fract   ! shift into first unit cell
REAL, DIMENSION(3) :: shift   ! shift into first unit cell
!
n_unit_cells(:) = cr_icc(:)
lout = .FALSE.
!
CALL do_readstru(strucfile, .FALSE.)
IF(ier_num/=0) THEN
   IF(ier_num /= 0) THEN
      IF(ier_msg(3) == ' ') THEN
         ier_msg(3) = strucfile(MAX(1,LEN_TRIM(strucfile)-LEN(ier_msg)):LEN_TRIM(strucfile))
      ENDIF 
   ENDIF 
   RETURN
ENDIF
CALL setup_lattice (cr_a0, cr_ar, cr_eps, cr_gten, cr_reps, &
            cr_rten, cr_win, cr_wrez, cr_v, cr_vr, lout, cr_gmat,       &
            cr_fmat, cr_cartesian,                                      &
            cr_tran_g, cr_tran_gi, cr_tran_f, cr_tran_fi)
CALL get_symmetry_matrices 
!
! Shift molecules such that the first atom has coordinates [0:1[
!
moles: DO im = 1, mole_num_mole
   vec(:)   = 0.0
   shift(:) = 0.0
   vec(:)   = cr_pos(:,mole_cont(mole_off(im) + 1))
   fract(1) = vec(1)   - REAL(INT(vec(1)))  + 1.0
   fract(1) = fract(1) - REAL(INT(fract(1)))
   shift(1) = fract(1) -vec(1)
   fract(2) = vec(2)   - REAL(INT(vec(2)))  + 1.0
   fract(2) = fract(2) - REAL(INT(fract(2)))
   shift(2) = fract(2) - vec(2)
   fract(3) = vec(3)   - REAL(INT(shift(3)))
   fract(3) = fract(3) - REAL(INT(fract(3)))
   shift(3) = fract(3) - vec(3)
   IF(MAXVAL(shift)>0.0) THEN
      atoms: DO j=1, mole_len(im)
         iat = mole_cont(mole_off(im) + j)
         cr_pos(1,iat) = cr_pos(1,iat) + shift(1)
         cr_pos(2,iat) = cr_pos(2,iat) + shift(2)
         cr_pos(3,iat) = cr_pos(3,iat) + shift(3)
      ENDDO atoms
   ENDIF
ENDDO moles
!
CALL save_store_setting             ! Backup user "save" setting
!ELLAcall save_show
CALL save_default_setting           ! Default to full saving
!ELLAcall save_show
line       = 'ignore, all'          ! Ignore all properties
length     = 11
CALL property_select(line, length, sav_sel_prop)
tempfile = 'internal'//'RBN_READMOLE'
CALL save_internal(tempfile)
CALL no_error
cr_icc(:) = n_unit_cells(:)         ! Restore intended number of unit cells
CALL readcell_internal(tempfile)
CALL save_restore_setting
CALL store_remove_single(tempfile, ier_num)
!
END SUBROUTINE readcell_mole
!
!*******************************************************************************
!
subroutine read_to_internal(infile, prefix)
!-
!   Reads the structure in 'infile' and stores it as prefix//infile in the 
!   internal storage
!   The current structure is lost !
!   Anny calling routine must preserver the old structure!
!+
!
use crystal_mod
use discus_save_mod
use prop_para_func
use prop_para_mod
use save_menu, ONLY: save_internal, save_store_setting, save_restore_setting, save_default_setting, save_struc, save_show
!
use precision_mod
!
implicit none
!
character(len=*), intent(in) :: infile
character(len=*), intent(in) :: prefix
!
character(len=PREC_STRING) :: line
character(len=PREC_STRING) :: getfile
integer :: length
logical :: l_site
!
call rese_cr
l_site = .false.
!
getfile = infile                    ! Just a local copy
call do_readstru(getfile, l_site)   ! Read actual file
if(ier_num/=0) return
!
call save_store_setting             ! Backup user "save" setting
if(ier_num/=0) return
!
call save_default_setting           ! Default to full saving
if(ier_num/=0) return
!
line       = 'ignore, all'          ! Ignore all properties
length     = 11
call property_select(line, length, sav_sel_prop)
if(ier_num/=0) return
!
line       = 'ignore, all'          ! Ignore all properties for global as well
length     = 11
call property_select(line, length,  cr_sel_prop)
if(ier_num/=0) return
!                                   ! Prepend with prefix
getfile = prefix(1:len_trim(prefix)) // infile(1:len_trim(infile))
!
call save_internal(getfile)         ! Finally save
if(ier_num/=0) return
!
call save_restore_setting           ! And restore user "save" settings
if(ier_num/=0) return
!
end subroutine read_to_internal
!
!*******************************************************************************
!
LOGICAL FUNCTION spcgr_test(spcgr)
!
USE spcgr_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN) :: spcgr
!
INTEGER        :: i
!
spcgr_test = .FALSE.
main: DO I=1, SPCGR_MAX
   IF(spcgr == spcgr_name(i)) THEN
      spcgr_test = .TRUE.
      EXIT main
   ENDIF
ENDDO main
!
END FUNCTION spcgr_test
!
!*******************************************************************************
!
END MODULE structur
