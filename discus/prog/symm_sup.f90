MODULE symm_sup_mod
!
CONTAINS
!*****7*****************************************************************
!
SUBROUTINE symm_setup 
!-                                                                      
!     Performs the generalized symmetry operation                       
!     See Sands, D.E. Vectors and Tensors in Crystallography Chapt. 4.7 
!+                                                                      
USE discus_config_mod 
use discus_allocate_appl_mod
USE crystal_mod 
USE metric_mod
use molecule_mod
USE symm_mod 
USE wyckoff_mod
!
USE errlist_mod 
USE trig_degree_mod
use precision_mod
!
IMPLICIT none 
!                                                                       
INTEGER :: i, j, k, l 
integer :: nscat
integer :: nsite
REAL(kind=PREC_DP) :: length 
!                                                                       
REAL(kind=PREC_DP) :: uij 
REAL(kind=PREC_DP) :: ctheta, stheta 
REAL(kind=PREC_DP) :: sym_d (3), sym_r (3) 
REAL(kind=PREC_DP) :: usym (4), ures (4) 
REAL(kind=PREC_DP) :: kron (3, 3) 
REAL(kind=PREC_DP) :: a (3, 3) 
REAL(kind=PREC_DP) :: b (3, 3) 
!                                                                       
!                                                                       
!                                                                       
DATA kron / 1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0 / 
!
nsite = 1
!
IF( cr_nscat > SYM_MAXSCAT .or. mole_num_type > SYM_MAXSCAT) THEN
   nscat = max ( cr_nscat, mole_num_type)
   nsite = max ( nsite, cr_ncatoms, SYM_MAXSITE)
   CALL alloc_symmetry ( nscat, nsite )
   IF ( ier_num < 0 ) THEN
      RETURN
   ENDIF
ENDIF
!write(*,*) 'SYM_UVW ', sym_uvw
!write(*,*) 'SYM_hkl ', sym_hkl
IF(sym_use == 0) THEN
!
!  Evaluate and use current EXPR
!
   call symm_check_expr
!                                                                       
!     initialize matrix and angle                                       
!
   sym_mat  = 0.0d0
   sym_rmat = 0.0d0
   sym_mat (4, 4) = 1.0d0 
   sym_rmat (4, 4) = 0.0d0 
!
   IF (sym_power_mult) then 
      ctheta = cosd (sym_angle) 
      stheta = sind (sym_angle) 
   ELSE 
      ctheta = cosd (sym_angle * sym_power) 
      stheta = sind (sym_angle * sym_power) 
   ENDIF 
!
!     symmetry origin is specified as absolute coordinates, or
!     via atom (from crystal or molecule)
!
   IF(sym_orig_type == 1) THEN       ! Origin defined by atom
      sym_orig(:) = cr_pos(:,sym_orig_atom)
   ENDIF
!
!     symmetry axis is specified as absolute coordinates, or 
!     via atom pair (from crystal or molecule)
!
   IF(sym_axis_type == 1) THEN       ! axis defined by atom pair
      DO i=1, 3
         sym_uvw(i) = cr_pos(i,sym_axis_atoms(2)) - cr_pos(i,sym_axis_atoms(1))
      ENDDO
      sym_hkl = matmul(real(cr_gten,KIND=PREC_DP), sym_uvw)
   ENDIF
!                                                                       
!     Create vectors of unit length in direct and reciprocal space      
!                                                                       
   length = sqrt(dot_product(sym_uvw, matmul(real(cr_gten, kind=PREC_DP), sym_uvw)))
   IF (length.eq.0.0) then 
      ier_num = - 32 
      ier_typ = ER_APPL 
      RETURN 
   ENDIF 
   DO j = 1, 3 
      sym_d (j) = sym_uvw (j) / length 
   ENDDO 
!
   length = sqrt(dot_product(sym_hkl, matmul(real(cr_rten, kind=PREC_DP), sym_hkl)))
   IF (length.eq.0.0) then 
      ier_num = - 32 
      ier_typ = ER_APPL 
      RETURN 
   ENDIF 
   DO j = 1, 3 
      sym_r (j) = sym_hkl (j) / length 
   ENDDO 
!                                                                       
!     calculate symmetry operation                                      
!                                                                       
   DO i = 1, 3 
      DO j = 1, 3 
         sym_mat (i, j) = 0.0d0 
         DO k = 1, 3 
            DO l = 1, 3 
               sym_mat(i, j) = sym_mat(i, j) + real(cr_rten(i, k), kind=PREC_DP) * &
                                               real(cr_eps(k, l,j),kind=PREC_DP) * sym_d(l)
            ENDDO 
         ENDDO 
         uij = sym_d (i) * sym_r (j) 
         sym_mat(i, j) = sym_mat(i, j) * stheta + uij + (kron(i, j) - uij) * ctheta
      ENDDO 
      IF (sym_power_mult) then 
         sym_mat(i, 4) = sym_trans (i) 
      ELSE 
         sym_mat(i, 4) = sym_trans (i) * sym_power 
      ENDIF 
   ENDDO 
!                                                                       
!     In case of improper rotation, multiply by -1                      
!                                                                       
   IF (.not.sym_type.and. (                                          &
      sym_power_mult.or..not.sym_power_mult.and.mod(sym_power, 2).ne.0) ) then
      DO i = 1, 3 
         DO j = 1, 3 
            sym_mat (i, j) = - sym_mat (i, j) 
         ENDDO 
      ENDDO 
   ENDIF 
!                                                                       
!     Transform symmetry operation into reciprocal space                
!                                                                       
   DO i = 1, 3 
      DO j = 1, 3 
         a (i, j) = sym_mat (i, j) 
      ENDDO 
   ENDDO 
!                                                                       
!     do transformation q = gSg*                                        
!                                                                       
   b = matmul(a, real(cr_rten, kind=PREC_DP))
   a = matmul(   real(cr_gten, kind=PREC_DP), b)
   DO i = 1, 3 
      DO j = 1, 3 
         sym_rmat (i, j) = a (i, j) 
      ENDDO 
   ENDDO 
!                                                                       
!     ----Calculate translational component due to origin               
!                                                                       
!                                                                       
!     ----Apply rotational part of symmetry operation to -1*origin      
!                                                                       
   DO j = 1, 3 
      usym (j) = - sym_orig (j) 
   ENDDO 
!                                                                       
!-----      ----Apply symmetry operation                                
!                                                                       
   usym (4) = 0.0 
!  CALL trans (usym, sym_mat, ures, 4) 
   ures = matmul(sym_mat,real(usym,KIND=PREC_DP))
!                                                                       
!     ----Add origin to result to obtain total translational part       
!                                                                       
   DO j = 1, 3 
      sym_or_tr (j) = ures (j) + sym_orig (j) 
   ENDDO 
ELSE
   sym_mat (:, :) = 0.0 
   sym_rmat(:, :) = 0.0 
   sym_mat (4, 4) = 1.0 
   sym_orig(:)    = 0.0
   sym_orig_mol   = .false.
!
   sym_mat (:, :) = spc_mat(:,:, sym_use)
!                                                                       
!     Transform symmetry operation into reciprocal space                
!                                                                       
   a(:,:) = sym_mat(1:3,1:3)
!                                                                       
!  do transformation q = gSg*                                        
!                                                                       
!  CALL matmulx (b, a, cr_rten) 
!  CALL matmulx (a, cr_gten, b) 
   b = matmul(a, real(cr_rten, kind=PREC_DP))
   a = matmul(   real(cr_gten, kind=PREC_DP), b)
   sym_rmat(1:3, 1:3) = a (:, :) 
!
ENDIF
!     write (output_io,2000) ((sym_rmat(i,j),j=1,3),i=1,3)              
!2000      format(3(3(2x,f10.6)/))                                      
!                                                                       
      END SUBROUTINE symm_setup                     
!
!*****7*****************************************************************
!
subroutine symm_check_expr
!-
!  Set the values of EXPR into the appropriate variables
!+
!
use symm_mod
!
use ber_params_mod
use do_replace_expr_mod
use errlist_mod
use precision_mod
!
integer, parameter :: MAXW = 13
integer :: ianz
character(len=PREC_STRING), dimension(MAXW) :: cpara
integer                   , dimension(MAXW) :: lpara
real(kind=PREC_DP)        , dimension(MAXW) :: werte
integer :: i
!
ianz = 13
do i=1,ianz
   cpara(i) = symm_expr(i)        ! Transfer into cpara, which may be modified
   lpara(i) = len_trim(cpara(i))
enddo
!
call do_use_expr(ianz, cpara, lpara, MAXW)        ! Use evaluated EXPR; EXPR string is replaced by numbers
call ber_params(ianz, cpara, lpara, werte, MAXW)  ! Evaluate, with current "EXPR"
if(ier_num/=0)  then
   ier_msg(1) = 'Expression for orient erroneous '
   return
endif
!
!   Replace sym* variable if EXPR was set. Test is not necessary, as defaults 
!   can be used if user might not have set the corresponding variable. 
if(any(symm_use_expr( 1: 3))) sym_uvw   = werte( 1:3)
if(any(symm_use_expr( 4: 6))) sym_hkl   = werte( 4:6)
if(    symm_use_expr( 7   ) ) sym_angle = werte( 7)
if(any(symm_use_expr( 8:10))) sym_trans = werte( 8:10)
if(any(symm_use_expr(11:13))) sym_orig  = werte(11:13)
!
end subroutine symm_check_expr
!
!*****7*****************************************************************
!
      SUBROUTINE symm_op_mult 
!-                                                                      
!     Performs the actual symmetry operation, multiple copy version     
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE atom_env_mod 
      USE modify_mod
      USE symm_mod 
      USE errlist_mod 
USE precision_mod
      IMPLICIT none 
!                                                                       
      INTEGER, PARAMETER :: MAXW = 6
!                                                                       
      CHARACTER(4) name 
!                                                                       
      INTEGER i, j, k, l 
      INTEGER i_start, i_end 
!                                                                       
      REAL(KIND=PREC_DP), dimension(4)    :: usym (4), ures (4) 
      REAL(KIND=PREC_DP), DIMENSION(MAXW) :: werte! (5) 
      REAL(KIND=PREC_DP), DIMENSION(4)    :: offset
!                                                                       
      DATA usym / 0.0D0, 0.0D0, 0.0D0, 1.0D0 / 
      DATA werte / 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 1.00D0 / 
!
      offset(4) = 0.0D0
!                                                                       
!     Set the appropriate starting end ending number for the atoms      
!                                                                       
      i_start = sym_start 
      i_end = sym_end 
      IF (sym_incl.eq.'all ') then 
         i_end = cr_natoms 
      ELSEIF (sym_incl.eq.'env ') then 
         i_end = atom_env (0) 
      ENDIF 
!                                                                       
!     Apply symmetry operation to all atoms within selected range       
!                                                                       
loop_atoms: DO l = i_start, i_end 
   i = l 
   IF (sym_incl.eq.'env ') i = atom_env (l) 
!                                                                       
!     --Select atom if:                                                 
!       type has been selected                                          
!                                                                       
      IF (sym_latom (cr_iscat (1,i) ) ) then 
        IF (sym_incl.eq.'env ') THEN 
           i = atom_env (l) 
           offset(1:3) = atom_pos(1:3,l) - cr_pos(1:3,i)
        ELSE
           offset(1:3) = 0.0
        ENDIF
!                                                                       
!     ----Subtract origin and possible offset due to periodic boundary
!                                                                       
         DO j = 1, 3 
         usym (j) = cr_pos (j, i) - sym_orig (j)  + offset(j)
         ENDDO 
!                                                                       
!-----      ----Apply symmetry operation                                
!                                                                       
         usym (4) = 1.0 
loop_pow:DO k = 1, sym_power 
!        CALL trans (usym, sym_mat, ures, 4) 
         ures = matmul(sym_mat, real(usym,KIND=PREC_DP))
!                                                                       
!     ----Add origin                                                    
!                                                                       
         DO j = 1, 3 
         werte (j + 1) = ures (j) + sym_orig (j) - offset(j)
         ENDDO 
         DO j = 1, 3 
         usym (j) = ures (j) 
         ENDDO 
         IF(sym_occup) THEN
            IF( symm_occupied(werte, sym_radius) ) CYCLE loop_pow
         ENDIF

!                                                                       
!     ----Insert copy of atom or replace original atom by its image     
!                                                                       
         IF (sym_mode) then 
            name = cr_at_lis (cr_iscat (1,i) ) 
            werte (5) = cr_dw (cr_iscat (1,i) ) 
            CALL do_ins_atom (name, MAXW, werte) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
         ELSE 
            DO j = 1, 3 
            cr_pos (j, i) = werte (j + 1) 
            ENDDO 
            IF (sym_incl.eq.'env ') THEN  ! Update the environment
               atom_pos(1:3,l) = cr_pos(1:3,i) + offset(1:3)
            ENDIF 
         ENDIF 
!        DO j = 1, 3 
!        usym (j) = ures (j) 
!        ENDDO 
      ENDDO  loop_pow
   ENDIF 
ENDDO  loop_atoms
!                                                                       
END SUBROUTINE symm_op_mult                   
!*****7*****************************************************************
      SUBROUTINE symm_op_single 
!-                                                                      
!     Performs the actual symmetry operation, single result version     
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE atom_env_mod 
      USE modify_mod
use prop_para_mod
      USE symm_mod 
      USE errlist_mod 
USE precision_mod
      IMPLICIT none 
!                                                                       
      INTEGER, PARAMETER :: MAXW = 6
!                                                                       
      CHARACTER(4) name 
      INTEGER i, j, l 
      INTEGER i_start, i_end 
      REAL(kind=PREC_DP), DIMENSION(4) :: usym (4), ures (4) 
      REAL(KIND=PREC_DP), DIMENSION(MAXW) :: werte !(5) 
      REAL(kind=PREC_DP), DIMENSION(4) :: offset
!                                                                       
      DATA usym / 0.0D0, 0.0D0, 0.0D0, 1.0D0 / 
      DATA werte / 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0 / 
!
      offset(4) = 0.0D0
!                                                                       
!     Set the appropriate starting end ending number for the atoms      
!                                                                       
      i_start = sym_start 
      i_end = sym_end 
      IF (sym_incl.eq.'all ') then 
         i_end = cr_natoms 
      ELSEIF (sym_incl.eq.'env ') then 
         i_end = atom_env (0) 
      ENDIF 
!                                                                       
!     Apply symmetry operation to all atoms within selected range       
!                                                                       
loop_atoms: DO l = i_start, i_end 
   i = l 
   IF (sym_incl.eq.'env ') i = atom_env (l) 
!                                                                       
!     --Select atom if:                                                 
!       type has been selected                                          
!                                                                       
      if(btest(cr_prop(i), PROP_TEMP)) cycle loop_atoms    ! 
        IF (sym_latom (cr_iscat (1,i) ) ) then 
          IF (sym_incl.eq.'env ') THEN 
             i = atom_env (l) 
             offset(1:3) = atom_pos(1:3,l) - cr_pos(1:3,i)
          ELSE
             offset(1:3) = 0.0
          ENDIF
!                                                                       
!     ----Subtract origin                                               
!                                                                       
         DO j = 1, 3 
         usym (j) = cr_pos (j, i) - sym_orig (j) + offset(j)
         ENDDO 
!                                                                       
!-----      ----Apply symmetry operation                                
!                                                                       
         usym (4) = 1.0 
!        CALL trans (usym, sym_mat, ures, 4) 
         ures = matmul(sym_mat, real(usym,KIND=PREC_DP))
!                                                                       
!     ----Add origin                                                    
!                                                                       
         DO j = 1, 3 
         werte (j + 1) = ures (j) + sym_orig (j) -offset(j)
         ENDDO 
!
!        If requested test for occupied positions
!
         IF(sym_occup) THEN
            IF( symm_occupied(werte, sym_radius) ) CYCLE loop_atoms
         ENDIF
!                                                                       
!     ----Insert copy of atom or replace original atom by its image     
!                                                                       
         IF (sym_mode) then 
            name = cr_at_lis (cr_iscat (1,i) ) 
            werte (5) = cr_dw (cr_iscat (1,i) ) 
            CALL do_ins_atom (name, MAXW, werte) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
         ELSE 
            DO j = 1, 3 
            cr_pos (j, i) = werte (j + 1) 
            ENDDO 
            IF (sym_incl.eq.'env ') THEN  ! Update the environment
               atom_pos(1:3,l) = cr_pos(1:3,i) + offset(1:3)
            ENDIF 
         ENDIF 
   ENDIF 
ENDDO loop_atoms
!                                                                       
END SUBROUTINE symm_op_single                 
!*****7*****************************************************************
      SUBROUTINE symm_mole_mult 
!-                                                                      
!     Performs the actual symmetry operation, multiple copy version     
!     Operates on molecules                                             
!+                                                                      
      USE discus_allocate_appl_mod
      USE discus_config_mod 
      USE crystal_mod 
      USE modify_mod
      USE molecule_mod 
      USE molecule_func_mod 
      USE spcgr_apply, ONLY: mole_insert_current
      USE symm_mod 
      USE errlist_mod 
USE precision_mod
      IMPLICIT none 
!                                                                       
      INTEGER, PARAMETER :: MAXW = 6
!                                                                       
      CHARACTER(4) name 
!                                                                       
      INTEGER i, j, k, l, ii , ll
      INTEGER i_start, i_end 
      INTEGER  :: at_start, at_end 
      INTEGER :: imole, imole_s, imole_t=1
      INTEGER  :: i1, i2
      INTEGER  :: n_gene   ! Number of molecule generators
      INTEGER  :: n_symm   ! Number of molecule symmetry operators
      INTEGER  :: n_mole   ! Number of molecules
      INTEGER  :: n_type   ! Number of molecule types
      INTEGER  :: n_atom   ! Number of atoms in molecules
      INTEGER  :: n_sub    ! Number of atoms in sub-molecules
      INTEGER  :: jatom    ! First atom in sub-molecule
      INTEGER, DIMENSION(:), ALLOCATABLE :: sub_list
      INTEGER, DIMENSION(:), ALLOCATABLE :: excl
!                                                                       
      REAL(KIND=PREC_DP) :: usym (4), ures (4), use_orig (3) 
      REAL(KIND=PREC_DP) , DIMENSION(MAXW) :: werte! (5) 
!                                                                       
      DATA usym / 0.0D0, 0.0D0, 0.0D0, 1.0D0 / 
      DATA werte / 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0 / 
!                                                                       
!     Set the appropriate starting end ending number for the atoms      
!                                                                       
      i_start = sym_start 
      i_end = sym_end 
      IF (sym_end.eq. - 1) i_end = mole_num_mole 
      imole_s = mole_num_mole 
!                                                                       
      IF (sym_new) then 
         n_gene = MAX( 1, MOLE_MAX_GENE)
         n_symm = MAX( 1, MOLE_MAX_SYMM)
         n_mole =         MOLE_MAX_MOLE
         n_type =         MOLE_MAX_TYPE
         n_atom =         MOLE_MAX_ATOM
         IF (mole_num_type+sym_power > MOLE_MAX_TYPE ) THEN
            n_type = MAX(mole_num_type + sym_power+5,MOLE_MAX_TYPE)
            call alloc_molecule(n_gene, n_symm, n_mole, n_type, n_atom)
         ENDIF
         imole_t = mole_num_type 
         IF (mole_num_type.lt.MOLE_MAX_TYPE) then 
            mole_num_type = mole_num_type+sym_power 
         ELSE 
            ier_num = - 66 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
      ENDIF 
!                                                                       
!     Apply symmetry operation to all molecules within selected range   
!                                                                       
      ll = 0
      DO l = i_start, i_end 
!                                                                       
!------ - Determine origin for symmetry operation                       
!                                                                       
      IF (sym_orig_mol) then 
         DO i = 1, 3 
         ii = mole_cont (mole_off (l) + 1) 
         use_orig (i) = sym_orig (i) + cr_pos (i, ii) 
         ENDDO 
      ELSE 
         DO i = 1, 3 
         use_orig (i) = sym_orig (i) 
         ENDDO 
      ENDIF 
!                                                                       
!     --Select atom if:                                                 
!       type has been selected                                          
!                                                                       
      IF (sym_latom (mole_type (l) ) ) then 
         ll = ll + 1
!                                                                       
!     ----Loop over all atoms in the molecule                           
!                                                                       
         IF(sym_sel_sub) THEN   ! Select a sub range of the molecule
            IF(ALLOCATED(sub_list)) DEALLOCATE(sub_list)
            ALLOCATE(sub_list(1:mole_len(l)))
            sub_list(:) = 0
            jatom = mole_cont (mole_off (l) +sym_sub_start)
            ALLOCATE(excl(1:sym_n_excl))
            DO i=1, sym_n_excl
               excl(i) = mole_cont(mole_off(l) + sym_excl(i))
            ENDDO
!write(*,*) ' MOLECULE NR ', l, ' jatom ', jatom,' Excl ',sym_n_excl,' : ', excl(:)
            CALL molecularize_sub(jatom,sym_n_excl,excl,mole_len(l),n_sub, sub_list)
            DEALLOCATE(excl)
            at_start = 1
            at_end   = n_sub
!write(*,*) ' MOLECULE NR ', l, ' Range ', at_start, at_end
         ELSE                   ! Select all atoms of molecule
            at_start = 1
            at_end   = mole_len(l)
         ENDIF
         IF(sym_axis_type == -1 .OR. sym_orig_type == -1) THEN  ! Need new setup
            IF(sym_axis_type == -1) THEN       ! axis defined by atom pair within molecule
               DO i=1, 3
                  i1 = mole_cont(mole_off(l)+sym_axis_atoms(1))
                  i2 = mole_cont(mole_off(l)+sym_axis_atoms(2))
                  sym_uvw(i) = cr_pos(i,i2) - cr_pos(i,i1)
               ENDDO
!              CALL trans (sym_uvw, cr_gten, sym_hkl, 3)
               sym_hkl = matmul(real(cr_gten,KIND=PREC_DP), sym_uvw)
            ENDIF
            IF(sym_orig_type == -1) THEN       ! origin defined by atom within molecule
               i1 = mole_cont(mole_off(l)+sym_orig_atom)
               sym_orig(:) = cr_pos(:,i1)
               use_orig(:) = sym_orig(:) 
            ENDIF
            CALL symm_setup
         ENDIF
         DO ii = at_start, at_end
            IF(sym_sel_sub) THEN   ! Select a sub range of the molecule
               i = sub_list(ii)
            ELSE
               i = mole_cont (mole_off (l) + ii) 
            ENDIF
!                                                                       
!     ----- Subtract origin                                             
!                                                                       
         DO j = 1, 3 
         usym (j) = cr_pos (j, i) - use_orig (j) 
         ENDDO 
!                                                                       
!-----      ----- Apply symmetry operation                              
!                                                                       
         usym (4) = 1.0 
         DO k = 1, sym_power 
!        CALL trans (usym, sym_mat, ures, 4) 
         ures = matmul(sym_mat, real(usym,KIND=PREC_DP))
!                                                                       
!     ------- Add origin                                                
!                                                                       
         DO j = 1, 3 
         werte (j + 1) = ures (j) + use_orig (j) 
         ENDDO 
!                                                                       
!     ------- Insert copy of atom or replace original atom by its image 
!                                                                       
         IF (sym_mode) then 
            name = cr_at_lis (cr_iscat (1,i) ) 
            werte (5) = cr_dw (cr_iscat (1,i) ) 
            CALL do_ins_atom (name, MAXW, werte) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
!                                                                       
!     --------- Insert atom into proper new molecule                    
!                                                                       
            imole = imole_s + (ll- i_start) * sym_power + (k - 1)       &
            + 1                                                         
            CALL mole_insert_current (cr_natoms, imole) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
         ELSE 
            DO j = 1, 3 
            cr_pos (j, i) = werte (j + 1) 
            ENDDO 
         ENDIF 
         DO j = 1, 3 
         usym (j) = ures (j) 
         ENDDO 
         ENDDO 
         ENDDO 
!                                                                       
!-----      ----Copy all properties                                     
!                                                                       
         IF (sym_mode) then 
            CALL copy_mole_char (mole_num_mole, l) 
         ENDIF 
!                                                                       
!----- ---- Set new molecule type if requested                          
!                                                                       
         IF (sym_new) then 
            DO k = 1, sym_power 
            IF (sym_mode) then 
               mole_type (mole_num_mole-sym_power + k) = imole_t + k 
            ELSE 
               mole_type (l) = imole_t + 1 
            ENDIF 
            ENDDO 
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
      END SUBROUTINE symm_mole_mult                 
!*****7*****************************************************************
      SUBROUTINE symm_mole_single 
!-                                                                      
!     Performs the actual symmetry operation, single result version     
!     Operates on molecules                                             
!+                                                                      
      USE discus_allocate_appl_mod
      USE discus_config_mod 
      USE crystal_mod 
      USE modify_mod
      USE molecule_mod 
      USE molecule_func_mod 
      USE spcgr_apply, ONLY: mole_insert_current
      USE symm_mod 
      USE errlist_mod 
USE precision_mod
      IMPLICIT none 
!                                                                       
      INTEGER, PARAMETER :: MAXW = 6
!                                                                       
      CHARACTER(4) name 
      INTEGER i, j, ii, l
      INTEGER i_start, i_end 
      INTEGER  :: at_start, at_end 
      INTEGER ::imole, imole_s, imole_t=1
      INTEGER  :: i1, i2
      INTEGER  :: n_gene   ! Number of molecule generators
      INTEGER  :: n_symm   ! Number of molecule symmetry operators
      INTEGER  :: n_mole   ! Number of molecules
      INTEGER  :: n_type   ! Number of molecule types
      INTEGER  :: n_atom   ! Number of atoms in molecules
      INTEGER  :: n_sub    ! Number of atoms in sub-molecules
      INTEGER  :: jatom    ! First atom in sub-molecule
      INTEGER, DIMENSION(:), ALLOCATABLE :: sub_list
      INTEGER, DIMENSION(:), ALLOCATABLE :: excl
      REAL(KIND=PREC_DP) :: usym (4), ures (4) 
      REAL(KIND=PREC_DP), DIMENSION(MAXW) :: werte !(5)
      REAL(kind=PREC_DP), DIMENSION(3)    :: use_orig !(3) 
!                                                                       
      DATA usym / 0.0D0, 0.0D0, 0.0D0, 1.0D0 / 
      DATA werte / 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0 , 1.0D0/ 
!                                                                       
!     Set the appropriate starting end ending number for the molecules  
!                                                                       
!write(*,*) ' MOLECULE SINGLE '
!write(*,*) ' start, end      ', sym_start, sym_end
!write(*,*) ' apply to type   ', sym_latom(1:2)
      i_start = sym_start 
      i_end = sym_end 
      IF (sym_end.eq. - 1) i_end = mole_num_mole 
      imole_s = mole_num_mole 
!                                                                       
      IF (sym_new) then 
         n_gene = MAX( 1, MOLE_MAX_GENE)
         n_symm = MAX( 1, MOLE_MAX_SYMM)
         n_mole =         MOLE_MAX_MOLE
         n_type =         MOLE_MAX_TYPE
         n_atom =         MOLE_MAX_ATOM
         IF (mole_num_type+sym_power > MOLE_MAX_TYPE ) THEN
            n_type = MAX(mole_num_type + sym_power+5,MOLE_MAX_TYPE)
            call alloc_molecule(n_gene, n_symm, n_mole, n_type, n_atom)
         ENDIF
         IF (mole_num_type.lt.MOLE_MAX_TYPE) then 
            mole_num_type = mole_num_type+1 
         ELSE 
            ier_num = - 66 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
         imole_t = mole_num_type 
      ENDIF 
!                                                                       
!     Apply symmetry operation to all molecules within selected range   
!                                                                       
      DO l = i_start, i_end 
!write(*,*) 'STARTING ON MOLE ', l, mole_num_mole
!                                                                       
!------ - Determine origin for symmetry operation                       
!                                                                       
      IF (sym_orig_mol) then 
         DO i = 1, 3 
         ii = mole_cont (mole_off (l) + 1) 
         use_orig (i) = sym_orig (i) + cr_pos (i, ii) 
         ENDDO 
      ELSE 
         DO i = 1, 3 
         use_orig (i) = sym_orig (i) 
         ENDDO 
      ENDIF 
!                                                                       
!     --Select atom if:                                                 
!       type has been selected                                          
!                                                                       
      IF (sym_latom (mole_type (l) ) ) then 
         IF (sym_mode) then 
            imole = mole_num_mole + 1
         ENDIF
!                                                                       
!     ----Loop over all atoms in the molecule                           
!                                                                       
         IF(sym_sel_sub) THEN   ! Select a sub range of the molecule
            IF(ALLOCATED(sub_list)) DEALLOCATE(sub_list)
            ALLOCATE(sub_list(1:mole_len(l)))
            sub_list(:) = 0
            jatom = mole_cont (mole_off (l) +sym_sub_start)
            ALLOCATE(excl(1:sym_n_excl))
            DO i=1, sym_n_excl
               excl(i) = mole_cont(mole_off(l) + sym_excl(i))
            ENDDO
            CALL molecularize_sub(jatom,sym_n_excl,excl,mole_len(l),n_sub, sub_list)
            DEALLOCATE(excl)
            at_start = 1
            at_end   = n_sub
         ELSE                   ! Select all atoms of molecule
            at_start = 1
            at_end   = mole_len(l)
         ENDIF
         IF(sym_axis_type == -1 .OR. sym_orig_type == -1) THEN  ! Need new setup
            IF(sym_axis_type == -1) THEN       ! axis defined by atom pair within molecule
               DO i=1, 3
                  i1 = mole_cont(mole_off(l)+sym_axis_atoms(1))
                  i2 = mole_cont(mole_off(l)+sym_axis_atoms(2))
                  sym_uvw(i) = cr_pos(i,i2) - cr_pos(i,i1)
               ENDDO
!              CALL trans (real(sym_uvw), cr_gten, real(sym_hkl), 3)
               sym_hkl = matmul(real(cr_gten,KIND=PREC_DP), sym_uvw)
            ENDIF
            IF(sym_orig_type == -1) THEN       ! origin defined by atom within molecule
               i1 = mole_cont(mole_off(l)+sym_orig_atom)
               sym_orig(:) = cr_pos(:,i1)
            ENDIF
            CALL symm_setup
         ENDIF
!write(*,*) 'Point A     MOLE ', l, mole_num_mole
         DO ii = at_start, at_end
            IF(sym_sel_sub) THEN   ! Select a sub range of the molecule
               i = sub_list(ii)
            ELSE
               i = mole_cont (mole_off (l) + ii) 
            ENDIF
!                                                                       
!     ----- Subtract origin                                             
!                                                                       
         DO j = 1, 3 
         usym (j) = cr_pos (j, i) - use_orig (j) 
         ENDDO 
!                                                                       
!-----      ----- Apply symmetry operation                              
!                                                                       
         usym (4) = 1.0 
!        CALL trans (usym, sym_mat, ures, 4) 
         ures = matmul(sym_mat,real(usym,KIND=PREC_DP))
!                                                                       
!     ----- Add origin                                                  
!                                                                       
         DO j = 1, 3 
         werte (j + 1) = ures (j) + use_orig (j) 
         ENDDO 
!                                                                       
!     ----- Insert copy of atom or replace original atom by its image   
!                                                                       
         IF (sym_mode) then 
            name = cr_at_lis (cr_iscat (1,i) ) 
            werte (5) = cr_dw (cr_iscat (1,i) ) 
            CALL do_ins_atom (name, MAXW, werte) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
!                                                                       
!     ------- Insert atom into proper new molecule                      
!                                                                       
!           imole = imole_s + l - i_start + 1 
!write(*,*) 'Point A 1   MOLE ', l, mole_num_mole, imole
            CALL mole_insert_current (cr_natoms, imole) 
!write(*,*) 'Point A 2   MOLE ', l, mole_num_mole, imole
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
         ELSE 
            DO j = 1, 3 
            cr_pos (j, i) = werte (j + 1) 
            ENDDO 
         ENDIF 
         ENDDO 
!write(*,*) 'Point B     MOLE ', l, mole_num_mole
!                                                                       
!-----      ----Copy all properties                                     
!                                                                       
         IF (sym_mode) then 
            CALL copy_mole_char (mole_num_mole, l) 
         ENDIF 
!                                                                       
!----- ---- Set new molecule type if requested                          
!                                                                       
         IF (sym_new) then 
            IF (sym_mode) then 
               mole_type (mole_num_mole) = imole_t 
            ELSE 
               mole_type (l) = imole_t 
            ENDIF 
         ENDIF 
      ENDIF 
!write(*,*) 'Point C     MOLE ', l, mole_num_mole
      ENDDO 
!write(*,*) 'FINISHED WITH ALL', mole_num_mole
!                                                                       
      END SUBROUTINE symm_mole_single               
!*****7*****************************************************************
      SUBROUTINE symm_domain_mult 
!-                                                                      
!     Performs the actual symmetry operation, multiple copy version     
!     Operates on molecules                                             
!+                                                                      
      USE discus_allocate_appl_mod
      USE discus_config_mod 
      USE crystal_mod 
      USE modify_mod
      USE molecule_mod 
      USE spcgr_apply, ONLY: mole_insert_current
      USE symm_mod 
      USE errlist_mod 
USE precision_mod
      IMPLICIT none 
!                                                                       
      INTEGER, PARAMETER :: MAXW = 6
!                                                                       
      CHARACTER(4) name 
!                                                                       
      INTEGER i, j, k, l, ii , ll
      INTEGER i_start, i_end 
      INTEGER :: imole, imole_s, imole_t=1
      INTEGER  :: n_gene   ! Number of molecule generators
      INTEGER  :: n_symm   ! Number of molecule symmetry operators
      INTEGER  :: n_mole   ! Number of molecules
      INTEGER  :: n_type   ! Number of molecule types
      INTEGER  :: n_atom   ! Number of atoms in molecules
!                                                                       
      REAL(KIND=PREC_DP) :: mat_atom (4, 4) 
      REAL(KIND=PREC_DP) :: mat_dime (4, 4) 
      REAL(KIND=PREC_DP) :: new_atom (4, 4) 
      REAL(KIND=PREC_DP) :: new_dime (4, 4) 
      REAL(KIND=PREC_DP) :: elements (3, 8) 
      REAL(KIND=PREC_DP) :: usym (4), ures (4), use_orig (3) 
      REAL(KIND=PREC_DP) , DIMENSION(MAXW) :: werte ! (5) 
!                                                                       
      DATA usym / 0.0, 0.0, 0.0, 1.0 / 
      DATA werte / 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0 / 
!                                                                       
!     Set the appropriate starting end ending number for the atoms      
!                                                                       
      i_start = sym_start 
      i_end = sym_end 
      IF (sym_end.eq. - 1) i_end = mole_num_mole 
      imole_s = mole_num_mole 
!                                                                       
      IF (sym_new) then 
         n_gene = MAX( 1, MOLE_MAX_GENE)
         n_symm = MAX( 1, MOLE_MAX_SYMM)
         n_mole =         MOLE_MAX_MOLE
         n_type =         MOLE_MAX_TYPE
         n_atom =         MOLE_MAX_ATOM
         IF (mole_num_type+sym_power > MOLE_MAX_TYPE ) THEN
            n_type = MAX(mole_num_type + sym_power+5,MOLE_MAX_TYPE)
            call alloc_molecule(n_gene, n_symm, n_mole, n_type, n_atom)
         ENDIF
         imole_t = mole_num_type 
         IF (mole_num_type.lt.MOLE_MAX_TYPE) then 
            mole_num_type = mole_num_type+sym_power 
         ELSE 
            ier_num = - 66 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
      ENDIF 
!                                                                       
!     Apply symmetry operation to all molecules within selected range   
!                                                                       
      ll = 0
      DO l = i_start, i_end 
!                                                                       
!------ - Determine origin for symmetry operation                       
!                                                                       
      IF (sym_orig_mol) then 
         DO i = 1, 3 
         ii = mole_cont (mole_off (l) + 1) 
         use_orig (i) = sym_orig (i) + cr_pos (i, ii) 
         ENDDO 
      ELSE 
         DO i = 1, 3 
         use_orig (i) = sym_orig (i) 
         ENDDO 
      ENDIF 
!                                                                       
!     --Select atom if:                                                 
!       type has been selected                                          
!                                                                       
      IF (sym_latom (mole_type (l) ) ) then 
         ll = ll + 1
!                                                                       
!     ----Create the matrices from the psueodoatom positions            
!                                                                       
         DO i = 1, 8 
         DO j = 1, 3 
         elements (j, i) = cr_pos (j, mole_cont (mole_off (l) + i) ) 
         ENDDO 
         ENDDO 
         DO i = 1, 3 
         DO j = 1, 3 
         mat_atom (i, j) = cr_pos (j, mole_cont (mole_off (l) + i + 1) )&
         - cr_pos (j, mole_cont (mole_off (l) + 1) )                    
         mat_dime (i, j) = cr_pos (j, mole_cont (mole_off (l) + i + 5) )&
         - cr_pos (j, mole_cont (mole_off (l) + 5) )                    
         ENDDO 
         mat_atom (4, i) = 0.0 
         mat_dime (4, i) = 0.0 
         mat_atom (i, 4) = cr_pos (i, mole_cont (mole_off (l) + 1) ) 
         mat_dime (i, 4) = cr_pos (i, mole_cont (mole_off (l) + 5) ) 
         ENDDO 
         mat_atom (4, 4) = 1.0 
         mat_dime (4, 4) = 1.0 
!                                                                       
!     ----Loop over Power of Operation                                  
!                                                                       
         DO k = 1, sym_power 
!        CALL matmul4 (new_atom, real(sym_mat), mat_atom) 
!        CALL matmul4 (new_dime, real(sym_mat), mat_dime) 
         new_atom = matmul(sym_mat, mat_atom)
         new_dime = matmul(sym_mat, mat_dime)
         DO i = 1, 3 
         DO j = 1, 3 
         mat_atom (i, j) = new_atom (i, j) 
         mat_dime (i, j) = new_dime (i, j) 
         ENDDO 
         ENDDO 
!                                                                       
!     ----Loop over all atoms in the molecule                           
!                                                                       
         DO ii = 1, 5, 4 
         i = mole_cont (mole_off (l) + ii) 
!                                                                       
!     ----- Subtract origin                                             
!                                                                       
         DO j = 1, 3 
         usym (j) = cr_pos (j, i) - use_orig (j) 
         usym (j) = elements (j, ii) - use_orig (j) 
         ENDDO 
!                                                                       
!-----      ----- Apply symmetry operation                              
!                                                                       
         usym (4) = 1.0 
!        CALL trans (usym, real(sym_mat), ures, 4) 
         ures = matmul(sym_mat,real(usym,KIND=PREC_DP))
!                                                                       
!     ------- Add origin                                                
!                                                                       
         DO j = 1, 3 
         werte (j + 1) = ures (j) + use_orig (j) 
         elements (j, ii) = ures (j) + use_orig (j) 
         ENDDO 
         ENDDO 
         IF (sym_dom_mode_atom) then 
            DO i = 1, 3 
            DO j = 1, 3 
            elements (j, 1 + i) = new_atom (i, j) + elements (j, 1) 
            ENDDO 
            ENDDO 
         ELSE 
            DO j = 1, 3 
            elements (j, 1) = cr_pos (j, mole_cont (mole_off (l)        &
            + 1) )                                                      
            ENDDO 
         ENDIF 
         IF (sym_dom_mode_shape) then 
            DO i = 1, 3 
            DO j = 1, 3 
            elements (j, 5 + i) = new_dime (i, j) + elements (j, 5) 
            ENDDO 
            ENDDO 
         ELSE 
            DO j = 1, 3 
            elements (j, 5) = cr_pos (j, mole_cont (mole_off (l)        &
            + 5) )                                                      
            ENDDO 
         ENDIF 
         DO ii = 1, 8 
         i = mole_cont (mole_off (l) + ii) 
         DO j = 1, 3 
         werte (j + 1) = elements (j, ii) 
         ENDDO 
!                                                                       
!     ------- Insert copy of atom or replace original atom by its image 
!                                                                       
         IF (sym_mode) then 
            name = cr_at_lis (cr_iscat (1,i) ) 
            werte (5) = cr_dw (cr_iscat (1,i) ) 
            CALL do_ins_atom (name, MAXW, werte) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
!                                                                       
!     --------- Insert atom into proper new molecule                    
!                                                                       
            imole = imole_s + (ll- i_start) * sym_power + (k - 1)       &
            + 1                                                         
            CALL mole_insert_current (cr_natoms, imole) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
         ELSE 
            DO j = 1, 3 
            cr_pos (j, i) = werte (j + 1) 
            ENDDO 
         ENDIF 
         DO j = 1, 3 
         usym (j) = ures (j) 
         ENDDO 
         ENDDO 
         ENDDO 
!                                                                       
!-----      ----Copy all properties                                     
!                                                                       
         IF (sym_mode) then 
            CALL copy_mole_char (mole_num_mole, l) 
         ENDIF 
!                                                                       
!----- ---- Set new molecule type if requested                          
!                                                                       
         IF (sym_new) then 
            DO k = 1, sym_power 
            IF (sym_mode) then 
               mole_type (mole_num_mole-sym_power + k) = imole_t + k 
            ELSE 
               mole_type (l) = imole_t + 1 
            ENDIF 
            ENDDO 
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
      END SUBROUTINE symm_domain_mult               
!*****7*****************************************************************
      SUBROUTINE symm_domain_single 
!-                                                                      
!     Performs the actual symmetry operation, single result version     
!     Operates on microdomain represenatation                           
!+                                                                      
      USE discus_allocate_appl_mod
      USE discus_config_mod 
      USE crystal_mod 
      USE modify_mod
      USE molecule_mod 
      USE spcgr_apply, ONLY: mole_insert_current
      USE symm_mod 
      USE errlist_mod 
USE precision_mod
      IMPLICIT none 
!                                                                       
      INTEGER, PARAMETER :: MAXW = 6
!                                                                       
      CHARACTER(len=4) name 
      INTEGER i, j, ii, l 
      INTEGER i_start, i_end 
      INTEGER :: imole, imole_s, imole_t=1
      INTEGER  :: n_gene   ! Number of molecule generators
      INTEGER  :: n_symm   ! Number of molecule symmetry operators
      INTEGER  :: n_mole   ! Number of molecules
      INTEGER  :: n_type   ! Number of molecule types
      INTEGER  :: n_atom   ! Number of atoms in molecules
      REAL(KIND=PREC_DP) :: mat_atom (4, 4) 
      REAL(KIND=PREC_DP) :: mat_dime (4, 4) 
      REAL(KIND=PREC_DP) :: new_atom (4, 4) 
      REAL(KIND=PREC_DP) :: new_dime (4, 4) 
      REAL(KIND=PREC_DP) :: elements (3, 8) 
      REAL(KIND=PREC_DP) :: usym (4), ures (4) 
      REAL(KIND=PREC_DP), DIMENSION(MAXW) :: werte! (5), use_orig (3) 
      REAL(KIND=PREC_DP), DIMENSION(3)    :: use_orig !(3) 
!                                                                       
      DATA usym / 0.0D0, 0.0D0, 0.0D0, 1.0D0 / 
      DATA werte / 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0 / 
!                                                                       
!     Set the appropriate starting end ending number for the molecules  
!                                                                       
      i_start = sym_start 
      i_end = sym_end 
      IF (sym_end.eq. - 1) i_end = mole_num_mole 
      imole_s = mole_num_mole 
!                                                                       
      IF (sym_new) then 
         n_gene = MAX( 1, MOLE_MAX_GENE)
         n_symm = MAX( 1, MOLE_MAX_SYMM)
         n_mole =         MOLE_MAX_MOLE
         n_type =         MOLE_MAX_TYPE
         n_atom =         MOLE_MAX_ATOM
         IF (mole_num_type+sym_power > MOLE_MAX_TYPE ) THEN
            n_type = MAX(mole_num_type + sym_power+5,MOLE_MAX_TYPE)
            call alloc_molecule(n_gene, n_symm, n_mole, n_type, n_atom)
         ENDIF
         IF (mole_num_type.lt.MOLE_MAX_TYPE) then 
            mole_num_type = mole_num_type+1 
         ELSE 
            ier_num = - 66 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
         imole_t = mole_num_type 
      ENDIF 
!                                                                       
!     Apply symmetry operation to all molecules within selected range   
!                                                                       
      DO l = i_start, i_end 
!                                                                       
!------ - Determine origin for symmetry operation                       
!                                                                       
      IF (sym_orig_mol) then 
         DO i = 1, 3 
         ii = mole_cont (mole_off (l) + 1) 
         use_orig (i) = sym_orig (i) + cr_pos (i, ii) 
         ENDDO 
      ELSE 
         DO i = 1, 3 
         use_orig (i) = sym_orig (i) 
         ENDDO 
      ENDIF 
!                                                                       
!     --Select atom if:                                                 
!       type has been selected                                          
!                                                                       
      IF (sym_latom (mole_type (l) ) ) then 
         IF(sym_mode) THEN
            imole = mole_num_mole + 1
         ENDIF
!                                                                       
!     ----Create the matrices from the psueodoatom positions            
!                                                                       
         DO i = 1, 8 
         DO j = 1, 3 
         elements (j, i) = cr_pos (j, mole_cont (mole_off (l) + i) ) 
         ENDDO 
         ENDDO 
         DO i = 1, 3 
         DO j = 1, 3 
         mat_atom (i, j) = cr_pos (j, mole_cont (mole_off (l) + i + 1) )&
         - cr_pos (j, mole_cont (mole_off (l) + 1) )                    
         mat_dime (i, j) = cr_pos (j, mole_cont (mole_off (l) + i + 5) )&
         - cr_pos (j, mole_cont (mole_off (l) + 5) )                    
         ENDDO 
         mat_atom (4, i) = 0.0 
         mat_dime (4, i) = 0.0 
         mat_atom (i, 4) = cr_pos (i, mole_cont (mole_off (l) + 1) ) 
         mat_dime (i, 4) = cr_pos (i, mole_cont (mole_off (l) + 5) ) 
         ENDDO 
         mat_atom (4, 4) = 1.0 
         mat_dime (4, 4) = 1.0 
!        CALL matmul4 (new_atom, real(sym_mat), mat_atom) 
!        CALL matmul4 (new_dime, real(sym_mat), mat_dime) 
         new_atom = matmul(sym_mat, mat_atom)
         new_dime = matmul(sym_mat, mat_dime)
!                                                                       
!     ----Loop over the two origins of the microdomain                  
!                                                                       
         DO ii = 1, 5, 4 
         i = mole_cont (mole_off (l) + ii) 
!                                                                       
!     ----- Subtract origin                                             
!                                                                       
         DO j = 1, 3 
         usym (j) = cr_pos (j, i) - use_orig (j) 
         ENDDO 
!                                                                       
!-----      ----- Apply symmetry operation                              
!                                                                       
         usym (4) = 1.0 
!        CALL trans (usym, real(sym_mat), ures, 4) 
         ures = matmul(sym_mat,real(usym,KIND=PREC_DP))
!                                                                       
!     ----- Add origin                                                  
!                                                                       
         DO j = 1, 3 
         werte (j + 1) = ures (j) + use_orig (j) 
         elements (j, ii) = ures (j) + use_orig (j) 
         ENDDO 
         ENDDO 
                                                                        
         IF (sym_dom_mode_atom) then 
            DO i = 1, 3 
            DO j = 1, 3 
            elements (j, 1 + i) = new_atom (i, j) + elements (j, 1) 
            ENDDO 
            ENDDO 
         ELSE 
            DO j = 1, 3 
            elements (j, 1) = cr_pos (j, mole_cont (mole_off (l)        &
            + 1) )                                                      
            ENDDO 
         ENDIF 
         IF (sym_dom_mode_shape) then 
            DO i = 1, 3 
            DO j = 1, 3 
            elements (j, 5 + i) = new_dime (i, j) + elements (j, 5) 
            ENDDO 
            ENDDO 
         ELSE 
            DO j = 1, 3 
            elements (j, 5) = cr_pos (j, mole_cont (mole_off (l)        &
            + 5) )                                                      
            ENDDO 
         ENDIF 
!                                                                       
!     ----- Insert copy of atom or replace original atom by its image   
!                                                                       
         DO ii = 1, 8 
         i = mole_cont (mole_off (l) + ii) 
         DO j = 1, 3 
         werte (j + 1) = elements (j, ii) 
         ENDDO 
         IF (sym_mode) then 
            name = cr_at_lis (cr_iscat (1,i) ) 
            werte (5) = cr_dw (cr_iscat (1,i) ) 
            CALL do_ins_atom (name, MAXW, werte) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
!                                                                       
!     ------- Insert atom into proper new molecule                      
!                                                                       
!           imole = imole_s + l - i_start + 1 
            CALL mole_insert_current (cr_natoms, imole) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
         ELSE 
            DO j = 1, 3 
            cr_pos (j, i) = werte (j + 1) 
            ENDDO 
         ENDIF 
         ENDDO 
!                                                                       
!-----      ----Copy all properties                                     
!                                                                       
         IF (sym_mode) then 
            CALL copy_mole_char (mole_num_mole, l) 
         ENDIF 
!                                                                       
!----- ---- Set new molecule type if requested                          
!                                                                       
         IF (sym_new) then 
            IF (sym_mode) then 
               mole_type (mole_num_mole) = imole_t 
            ELSE 
               mole_type (l) = imole_t 
            ENDIF 
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
      END SUBROUTINE symm_domain_single             
!*****7*****************************************************************
!
SUBROUTINE symm_ca_mult(uvw, lspace, loutput) 
!-                                                                      
!     Performs the actual symmetry operation, multiple copy version     
!     Only the input vector uvw is used in direct or reciprocal space   
!+                                                                      
USE discus_config_mod 
USE symm_mod 
!                                                                       
USE errlist_mod 
USE param_mod 
USE precision_mod
USE prompt_mod 
!
IMPLICIT none 
!                                                                       
REAL(kind=PREC_DP)   , DIMENSION(3), INTENT(INOUT) :: uvw
LOGICAL,               INTENT(IN)    :: lspace 
LOGICAL,               INTENT(IN)    :: loutput 
!
INTEGER :: j, k 
!                                                                       
REAL(KIND=PREC_DP), DIMENSION(4) :: usym, ures 
REAL(KIND=PREC_DP), DIMENSION(5) :: werte
!                                                                       
DATA usym / 0.0D0, 0.0D0, 0.0D0, 1.0D0 / 
DATA werte / 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0 / 
!                                                                       
!     real space part                                                   
!                                                                       
IF (lspace) THEN 
!                                                                       
!     ----Subtract origin, if in real space                             
!                                                                       
   DO j = 1, 3 
      usym(j) = uvw(j) - sym_orig(j) 
   ENDDO 
!                                                                       
!-----      --Apply symmetry operation                                  
!                                                                       
   usym(4) = 1.0 
   DO k = 1, sym_power 
!     CALL trans (usym, real(sym_mat), ures, 4) 
      ures = matmul(sym_mat,real(usym,KIND=PREC_DP))
!                                                                       
!     ----Add origin and store result                                   
!                                                                       
      DO j = 1, 3 
         res_para((k - 1) * 3 + j) = ures(j) + sym_orig(j) 
      ENDDO 
      IF(loutput) THEN
         WRITE (output_io, 3000) (res_para( (k - 1) * 3 + j), j = 1, 3) 
      ENDIF
!                                                                       
!     ----Replace current vector by its image                           
!                                                                       
      DO j = 1, 3 
         usym(j) = ures(j) 
      ENDDO 
   ENDDO 
ELSE 
!                                                                       
!     ----Reciprocal space part                                         
!                                                                       
   DO j = 1, 3 
      usym(j) = uvw(j) 
   ENDDO 
!                                                                       
!-----      --Apply symmetry operation                                  
!                                                                       
   usym(4) = 0.0 
   DO k = 1, sym_power 
!     CALL trans (usym, real(sym_rmat), ures, 4) 
      ures = matmul(sym_rmat,real(usym,KIND=PREC_DP))
!                                                                       
!     ----Store result                                                  
!                                                                       
      DO j = 1, 3 
         res_para((k - 1) * 3 + j) = ures(j) 
      ENDDO 
      IF(loutput) THEN
         WRITE(output_io, 3000) (res_para((k - 1) * 3 + j), j = 1, 3) 
      ENDIF
!                                                                       
!     ----Replace current vector by its image                           
!                                                                       
      DO j = 1, 3 
         usym(j) = ures(j) 
      ENDDO 
   ENDDO 
ENDIF 
!                                                                       
res_para (0) = sym_power * 3 
!                                                                       
 3000 FORMAT    (' Result    : ',3(2x,f9.4)) 
!
END SUBROUTINE symm_ca_mult                   
!*****7*****************************************************************
      SUBROUTINE symm_ca_single (uvw, lspace, loutput) 
!-                                                                      
!     Performs the actual symmetry operation, multiple copy version     
!     Only the input vector uvw is used in direct or reciprocal space   
!+                                                                      
      USE discus_config_mod 
      USE symm_mod 
!                                                                       
      USE errlist_mod 
      USE param_mod 
USE precision_mod
      USE prompt_mod 
!
      IMPLICIT none 
!                                                                       
      REAL(kind=PREC_DP)   ,DIMENSION(1:3), INTENT(IN) :: uvw
      LOGICAL,                INTENT(IN) :: lspace 
      LOGICAL,                INTENT(IN) :: loutput 
!                                                                       
      INTEGER j 
!
      REAL(kind=PREC_DP) :: usym (4), ures (4) 
      REAL(KIND=PREC_DP) :: werte (5) 
!                                                                       
      DATA usym / 0.0D0, 0.0D0, 0.0D0, 1.0D0 / 
      DATA werte / 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0 / 
!                                                                       
!     real space part                                                   
!                                                                       
      IF (lspace) then 
!                                                                       
!     ----Subtract origin, if in real space                             
!                                                                       
         DO j = 1, 3 
         usym (j) = uvw (j) - sym_orig (j) 
         ENDDO 
!                                                                       
!-----      --Apply symmetry operation                                  
!                                                                       
         usym (4) = 1.0 
!        CALL trans (usym, REAL(sym_mat), ures, 4) 
         ures = matmul(sym_mat,real(usym,KIND=PREC_DP))
!                                                                       
!     ----Add origin and store result                                   
!                                                                       
         DO j = 1, 3 
         res_para (j) = ures (j) + sym_orig (j) 
         ENDDO 
!                                                                       
!     ----Replace current vector by its image                           
!                                                                       
         DO j = 1, 3 
         usym (j) = ures (j) 
         ENDDO 
      ELSE 
!                                                                       
!     ----Subtract origin, if in real space                             
!                                                                       
         DO j = 1, 3 
         usym (j) = uvw (j) 
         ENDDO 
!                                                                       
!-----      --Apply symmetry operation                                  
!                                                                       
         usym (4) = 0.0 
!        CALL trans (usym, real(sym_rmat), ures, 4) 
         ures = matmul(sym_rmat,real(usym,KIND=PREC_DP))
!                                                                       
!     ----Add origin and store result                                   
!                                                                       
         DO j = 1, 3 
         res_para (j) = ures (j) 
         ENDDO 
!                                                                       
!     ----Replace current vector by its image                           
!                                                                       
         DO j = 1, 3 
         usym (j) = ures (j) 
         ENDDO 
      ENDIF 
!                                                                       
      IF( loutput) THEN
         res_para (0) = 3 
         WRITE (output_io, 3000) (res_para (j), j = 1, 3) 
      ENDIF
!                                                                       
 3000 FORMAT    (' Result    : ',3(2x,f9.4)) 
      END SUBROUTINE symm_ca_single                 
!
!*******************************************************************************
!
subroutine symm_rec
!-
! Perform the symmetry operations on the recursive environment of a given atom
! The 'include', 'select' status are temporarily stored. The lower routines 
! must include all atoms and all atom types that are within the current 
! environment.
! The temporary property flag is used to perform the operation only once, i.e.
! only on those atoms whose property is not set. At the end the property is
! cleared globally.
!+
!
use crystal_mod
use atom_env_mod
use do_find_mod
use prop_para_mod
use symm_mod
!
implicit none
!
integer, parameter :: MAXW = 1
integer :: i, j     ! Dummy index
integer :: ianz     ! Number of params in werte
integer :: r_start  ! Main symmetry inlude start
integer :: r_end    ! Main symmetry inlude end
!logical,          dimension(:), allocatable  ::  r_sym_latom   ! main symmetry selected atom list
!logical,          dimension(:), allocatable  ::  r_sym_lsite   ! main symmetry site list
real(kind=PREC_DP), dimension(MAXW) :: werte                   ! Array to transfer atom number to conn
!
sym_incl = 'env'                                   ! For lower routines use 'env' mode
r_start = sym_start                                ! Retain original include limits (start / end)
if(sym_end==-1) then
   r_end = cr_natoms                               ! User specified "all" on 'include' command
else
   r_end = sym_end                                 ! User specified limits
endif
!
do i=1, cr_natoms
   cr_prop(i) = ibclr(cr_prop(i), PROP_TEMP)       ! Set   temporary property for all atoms
enddo
cr_sel_prop(0) = ibclr(cr_sel_prop(0), PROP_TEMP)  ! Set temporary property as 'present'
cr_sel_prop(1) = ibclr(cr_sel_prop(1), PROP_TEMP)  ! Set temporary property as 'present'
!
loop_incl: do i=r_start, r_end                     ! Loop over all atoms included
   if(.not. sym_latom(cr_iscat(1,i))) cycle loop_incl      ! Only for selected atoms
   if(btest(cr_prop(i), PROP_TEMP)) cycle loop_incl
   ianz = 1
   werte(1) = real(i)
!  werte(2) =  1
   call do_find_neig_conn_all(ianz, MAXW, werte)
!write(*,*) ' CENTRAL ', i, btest(cr_prop(i), PROP_TEMP)
!write(*,*) ' neig    ', atom_env(1:atom_env(0))
   sym_start = 1
   sym_end   = atom_env(0)
!
!  Since each connectivity might/should have its own random values set up symmetry
!
   call symm_setup 
!
!  Now perform actual symmetry operation
!
   if(sym_power_mult) THEN
      call symm_op_mult
   else
      call symm_op_single
   endif
   do j=1, atom_env(0)
      cr_prop(atom_env(j)) = ibset(cr_prop(atom_env(j)), PROP_TEMP)  ! Set TEMP prop
   enddo
enddo loop_incl
!
do i=1, cr_natoms
   cr_prop(i) = ibclr(cr_prop(i), PROP_TEMP)       ! Clear temporary property
enddo
!
!
end subroutine symm_rec
!
!*******************************************************************************
!
LOGICAL FUNCTION symm_occupied(werte, radius)
!
USE crystal_mod
USE metric_mod
USE precision_mod
!
IMPLICIT NONE
!
REAL(KIND=PREC_DP), DIMENSION(5), INTENT(IN) :: werte
REAL(KIND=PREC_DP), INTENT(IN) :: radius
!
LOGICAL, PARAMETER    ::LSPACE = .TRUE.
INTEGER               :: i
REAL(kind=PREC_DP)   , DIMENSION(3) :: pos
REAL(kind=PREC_DP)   , DIMENSION(3) :: vec
!
pos(:) = werte(2:4)
symm_occupied = .FALSE.
main:DO i=1, cr_natoms
   vec(:) = cr_pos(:,i)
   symm_occupied = do_blen(LSPACE, pos,vec) <= radius
   IF(symm_occupied) EXIT main
ENDDO main
!
END FUNCTION symm_occupied
!
!*****7*****************************************************************
!
SUBROUTINE symm_store
!
USE discus_allocate_appl_mod
USE symm_mod
USE symm_temp_mod
!
IMPLICIT NONE
!
INTEGER :: i,j
!
SYM_TEMP_MAXSCAT = SYM_MAXSCAT
SYM_TEMP_MAXSITE = SYM_MAXSITE
!
IF(ALLOCATED(sym_latom)) THEN
   IF(ALLOCATED(sym_temp_latom)) DEALLOCATE(sym_temp_latom)
   i = LBOUND(sym_latom,1)
   j = UBOUND(sym_latom,1)
   ALLOCATE(sym_temp_latom(i:j))
   sym_temp_latom(i:j) = sym_latom(i:j)
ENDIF
!
IF(ALLOCATED(sym_lsite)) THEN
   IF(ALLOCATED(sym_temp_lsite)) DEALLOCATE(sym_temp_lsite)
   i = LBOUND(sym_lsite,1)
   j = UBOUND(sym_lsite,1)
   ALLOCATE(sym_temp_lsite(i:j))
   sym_temp_lsite(i:j) = sym_lsite(i:j)
ENDIF
!
IF(ALLOCATED(sym_excl)) THEN 
   IF(ALLOCATED(sym_temp_excl)) DEALLOCATE(sym_temp_excl)
   i = LBOUND(sym_excl,1)
   j = UBOUND(sym_excl,1)
   ALLOCATE(sym_temp_excl(i:j))
   sym_temp_excl(i:j) = sym_excl(i:j)
ENDIF
!
sym_temp_incl       = sym_incl 
sym_temp_use        = sym_use  ! ace group symmetry operation no. N
sym_temp_sel_mode   = sym_sel_mode
sym_temp_sel_prop   = sym_sel_prop
sym_temp_start      = sym_start
sym_temp_end        = sym_end
sym_temp_sub_start  = sym_sub_start
sym_temp_sub_end    = sym_sub_end
sym_temp_n_excl     = sym_n_excl
sym_temp_power      = sym_power
sym_temp_axis_type  = sym_axis_type
sym_temp_axis_atoms = sym_axis_atoms
sym_temp_orig_type  = sym_orig_type
sym_temp_orig_atom  = sym_orig_atom
sym_temp_mode       = sym_mode
sym_temp_new        = sym_new
sym_temp_orig_mol   = sym_orig_mol
sym_temp_power_mult = sym_power_mult
sym_temp_type       = sym_type
sym_temp_occup      = sym_occup
sym_temp_sel_atom   = sym_sel_atom
sym_temp_sel_sub    = sym_sel_sub
sym_temp_dom_mode_shape = sym_dom_mode_shape
sym_temp_dom_mode_atom  = sym_dom_mode_atom
sym_temp_angle     = sym_angle
sym_temp_radius    = sym_radius 
sym_temp_hkl       = sym_hkl
sym_temp_orig      = sym_orig
sym_temp_or_tr     = sym_or_tr
sym_temp_trans     = sym_trans
sym_temp_uvw       = sym_uvw
sym_temp_mat       = sym_mat
sym_temp_rmat      = sym_rmat
!
END SUBROUTINE symm_store
!
!*****7*****************************************************************
!
SUBROUTINE symm_restore
!
USE discus_allocate_appl_mod
USE symm_mod
USE symm_temp_mod
!
IMPLICIT NONE
!
INTEGER :: i,j
!
SYM_MAXSCAT = SYM_TEMP_MAXSCAT
SYM_MAXSITE = SYM_TEMP_MAXSITE
!
IF(ALLOCATED(sym_temp_latom)) THEN
   IF(ALLOCATED(sym_latom)) DEALLOCATE(sym_latom)
   i = LBOUND(sym_temp_latom,1)
   j = UBOUND(sym_temp_latom,1)
   ALLOCATE(sym_latom(i:j))
   sym_latom(i:j) = sym_temp_latom(i:j)
ENDIF
!
IF(ALLOCATED(sym_temp_lsite)) THEN
   IF(ALLOCATED(sym_lsite)) DEALLOCATE(sym_lsite)
   i = LBOUND(sym_temp_lsite,1)
   j = UBOUND(sym_temp_lsite,1)
   ALLOCATE(sym_lsite(i:j))
   sym_lsite(i:j) = sym_temp_lsite(i:j)
ENDIF
!
IF(ALLOCATED(sym_temp_excl)) THEN
   IF(ALLOCATED(sym_excl)) DEALLOCATE(sym_excl)
   i = LBOUND(sym_temp_excl,1)
   j = UBOUND(sym_temp_excl,1)
   ALLOCATE(sym_excl(i:j))
   sym_excl(i:j) = sym_temp_excl(i:j)
ENDIF
!
sym_incl       = sym_temp_incl 
sym_use        = sym_temp_use  ! ace group symmetry operation no. N
sym_sel_mode   = sym_temp_sel_mode
sym_sel_prop   = sym_temp_sel_prop
sym_start      = sym_temp_start
sym_end        = sym_temp_end
sym_sub_start  = sym_temp_sub_start
sym_sub_end    = sym_temp_sub_end
sym_n_excl     = sym_temp_n_excl
sym_power      = sym_temp_power
sym_axis_type  = sym_temp_axis_type
sym_axis_atoms = sym_temp_axis_atoms
sym_orig_type  = sym_temp_orig_type
sym_orig_atom  = sym_temp_orig_atom
sym_mode       = sym_temp_mode
sym_new        = sym_temp_new
sym_orig_mol   = sym_temp_orig_mol
sym_power_mult = sym_temp_power_mult
sym_type       = sym_temp_type
sym_occup      = sym_temp_occup
sym_sel_atom   = sym_temp_sel_atom
sym_sel_sub    = sym_temp_sel_sub
sym_dom_mode_shape = sym_temp_dom_mode_shape
sym_dom_mode_atom  = sym_temp_dom_mode_atom
sym_angle     = sym_temp_angle
sym_radius    = sym_temp_radius 
sym_hkl       = sym_temp_hkl
sym_orig      = sym_temp_orig
sym_or_tr     = sym_temp_or_tr
sym_trans     = sym_temp_trans
sym_uvw       = sym_temp_uvw
sym_mat       = sym_temp_mat
sym_rmat      = sym_temp_rmat
!
END SUBROUTINE symm_restore
!
!*****7*****************************************************************
!
SUBROUTINE symm_reset
!
USE discus_allocate_appl_mod
USE symm_mod
!
IMPLICIT NONE
!
CALL alloc_symmetry(1, 1)
!
IF(ALLOCATED(sym_latom)) sym_latom(:) = .FALSE.   ! (0:SYM_MAXSCAT)
IF(ALLOCATED(sym_lsite)) sym_lsite(:) = .TRUE.    ! (0:SYM_MAXSCAT)
IF(ALLOCATED(sym_excl))  sym_excl     = 0
sym_incl       = 'list'
sym_use        = 0  ! Use space group symmetry operation no. N
sym_sel_mode   = 0
sym_sel_prop   = (/0,0/)
sym_start      = 1
sym_end        = 1
sym_sub_start  = 1
sym_sub_end    = 1
sym_n_excl     = 1
sym_power      = 1
sym_axis_type  = 0  ! axis type (0 abolute, 1 atoms in crystal, -1 atoms in mol
sym_axis_atoms = 0  ! Atoms that define the axis
sym_orig_type  = 0  ! origin type (0 abolute, 1 atoms in crystal, -1 atoms in mol
sym_orig_atom  = 1  ! Atom at origin of symmetry operation
sym_mode       = .true.
sym_new        = .false.
sym_orig_mol   = .true.
sym_power_mult = .true.
sym_type       = .true.
sym_occup      = .false.
sym_occup      = .false.
sym_sel_atom   = .true.
sym_sel_sub    = .FALSE.
sym_dom_mode_shape = .false.
sym_dom_mode_atom  = .false.
sym_angle     = 0.0
sym_radius    = 1.0E-8
sym_hkl       = (/0.0, 0.0, 1.0/)
sym_orig      = (/0.0, 0.0, 0.0/)
sym_or_tr     = (/0.0, 0.0, 0.0/)
sym_trans     = (/0.0, 0.0, 0.0/)
sym_uvw       = (/0.0, 0.0, 1.0/)
sym_mat       = 0.0
sym_rmat      = 0.0
!
END SUBROUTINE symm_reset
!
!*****7*****************************************************************
!
END MODULE symm_sup_mod
