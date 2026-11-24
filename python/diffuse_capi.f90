MODULE diffuse_capi
!
!  C interface to the DISCUS suite (migration plan WP2, option B; see
!  migration/WP2_Binding_Decision.md and python/diffuse_core.h).
!
!  Plain C types only, so that any language can call it (Python through
!  ctypes, C, Julia ...):
!    * strings in : char array + 64-bit length (no terminating zero needed)
!    * strings out: char buffer of a fixed, documented size, blank padded
!    * arrays     : pointer + 64-bit size, explicit-size in Fortran
!    * sizes      : 64-bit; flags and codes: int
!  The C names start with "dc_"; the Fortran names here start with "capi_"
!  (the module "suite" in diffuse_suite.f90 already has dc_init etc.).
!
!  Everything runs in library mode (dc_init): an error never ends the
!  process; the first error is kept for dc_last_error.
!
USE, INTRINSIC :: ISO_C_BINDING
USE suite
!
IMPLICIT NONE
!
PRIVATE
!
INTEGER(C_INT), PARAMETER :: API_VERSION = 4   ! 2: crystal out, 3: in, 4: ADPs, molecules (WP3)
!
CONTAINS
!
!*******************************************************************************
!
FUNCTION from_c(text, n) RESULT(s)
!
!  A Fortran string from a C char array of length n
!
INTEGER(C_INT64_T),                      INTENT(IN) :: n
CHARACTER(KIND=C_CHAR), DIMENSION(n),    INTENT(IN) :: text
CHARACTER(LEN=n) :: s
INTEGER(C_INT64_T) :: i
!
DO i = 1, n
   s(i:i) = text(i)
ENDDO
!
END FUNCTION from_c
!
!*******************************************************************************
!
SUBROUTINE to_c(s, buf, n)
!
!  Copy a Fortran string into a C buffer of size n, blank padded
!
CHARACTER(LEN=*),                        INTENT(IN)  :: s
INTEGER(C_INT64_T),                      INTENT(IN)  :: n
CHARACTER(KIND=C_CHAR), DIMENSION(n),    INTENT(OUT) :: buf
INTEGER(C_INT64_T) :: i
!
buf = ' '
DO i = 1, MIN(n, INT(LEN(s), C_INT64_T))
   buf(i) = s(i:i)
ENDDO
!
END SUBROUTINE to_c
!
!*******************************************************************************
!
FUNCTION capi_api_version() BIND(C, NAME='dc_api_version') RESULT(v)
INTEGER(C_INT) :: v
v = API_VERSION
END FUNCTION capi_api_version
!
!*******************************************************************************
!
SUBROUTINE capi_init() BIND(C, NAME='dc_init')
CALL dc_init()
END SUBROUTINE capi_init
!
!*******************************************************************************
!
SUBROUTINE capi_reset() BIND(C, NAME='dc_reset')
CALL dc_reset()
END SUBROUTINE capi_reset
!
!*******************************************************************************
!
SUBROUTINE capi_flush_output() BIND(C, NAME='dc_flush_output')
CALL flush_output()
END SUBROUTINE capi_flush_output
!
!*******************************************************************************
!
FUNCTION capi_run(text, n) BIND(C, NAME='dc_run') RESULT(longest)
!
!  Run command lines (separated by newlines) as a macro held in memory.
!  Returns 0 when they ran; otherwise (nothing run) the length of the
!  longest line, which exceeds the limit of 1024 characters.
!
INTEGER(C_INT64_T),                      VALUE      :: n
CHARACTER(KIND=C_CHAR), DIMENSION(n),    INTENT(IN) :: text
INTEGER(C_INT64_T) :: longest
INTEGER :: nchar_max, nlines
!
IF(n == 0) THEN
   longest = 0
   RETURN
ENDIF
CALL dc_run_commands(from_c(text, n), nchar_max, nlines)
longest = 0
IF(nchar_max > 1024) longest = nchar_max
!
END FUNCTION capi_run
!
!*******************************************************************************
!
FUNCTION capi_value(expr, n) BIND(C, NAME='dc_value') RESULT(v)
!
!  Value of a DISCUS expression, as a double (about 12 significant digits,
!  as many as the expression evaluator keeps). On an error the result is 0
!  and dc_last_error tells what went wrong (clear it first).
!
USE ber_params_mod
USE errlist_mod
USE lib_errlist_func
USE precision_mod
!
INTEGER(C_INT64_T),                      VALUE      :: n
CHARACTER(KIND=C_CHAR), DIMENSION(n),    INTENT(IN) :: expr
REAL(C_DOUBLE) :: v
CHARACTER(LEN=PREC_STRING), DIMENSION(1) :: cpara
INTEGER                   , DIMENSION(1) :: lpara
REAL(KIND=PREC_DP)        , DIMENSION(1) :: werte
!
v = 0.0_C_DOUBLE
IF(n < 1 .OR. n > PREC_STRING) RETURN
cpara(1) = from_c(expr, n)
lpara(1) = LEN_TRIM(cpara(1))
werte    = 0.0_PREC_DP
CALL ber_params(1, cpara, lpara, werte, 1)
IF(ier_num /= 0) THEN
   CALL errlist
   RETURN
ENDIF
v = REAL(werte(1), C_DOUBLE)
!
END FUNCTION capi_value
!
!*******************************************************************************
!
SUBROUTINE capi_set_value(var, n, val) BIND(C, NAME='dc_set_value')
!
!  Set a user variable (declared with "variable real, name"). Built-in
!  arrays such as r[] are set with a command (dc_run "r[2] = ...").
!
USE do_variable_mod
USE errlist_mod
USE lib_errlist_func
USE precision_mod
!
INTEGER(C_INT64_T),                      VALUE      :: n
CHARACTER(KIND=C_CHAR), DIMENSION(n),    INTENT(IN) :: var
REAL(C_DOUBLE),                          VALUE      :: val
CHARACTER(LEN=PREC_STRING) :: name
CHARACTER(LEN=1)           :: dummy
REAL(KIND=PREC_DP)         :: wert
INTEGER                    :: laenge, length
INTEGER, DIMENSION(2)      :: substr
!
IF(n < 1 .OR. n > PREC_STRING) RETURN
name   = from_c(var, n)
laenge = LEN_TRIM(name)
dummy  = ' '
wert   = REAL(val, PREC_DP)
CALL upd_variable(name, laenge, wert, dummy, length, substr)
IF(ier_num /= 0) CALL errlist
!
END SUBROUTINE capi_set_value
!
!*******************************************************************************
!
SUBROUTINE capi_clear_error() BIND(C, NAME='dc_clear_error')
CALL dc_clear_error()
END SUBROUTINE capi_clear_error
!
!*******************************************************************************
!
SUBROUTINE capi_last_error(num, typ, cls, text, msgs) BIND(C, NAME='dc_last_error')
!
!  The first error since dc_clear_error: number (< 0; 0 if none), type,
!  type name (4 characters), standard message (45) and DISCUS's extra
!  lines (7 x 80 = 560 characters), all blank padded.
!
INTEGER(C_INT),                           INTENT(OUT) :: num
INTEGER(C_INT),                           INTENT(OUT) :: typ
CHARACTER(KIND=C_CHAR), DIMENSION(4),     INTENT(OUT) :: cls
CHARACTER(KIND=C_CHAR), DIMENSION(45),    INTENT(OUT) :: text
CHARACTER(KIND=C_CHAR), DIMENSION(560),   INTENT(OUT) :: msgs
INTEGER            :: fnum, ftyp
CHARACTER(LEN=4)   :: fcls
CHARACTER(LEN=45)  :: ftext
CHARACTER(LEN=560) :: fmsgs
!
CALL dc_last_error(fnum, ftyp, fcls, ftext, fmsgs)
num = fnum
typ = ftyp
CALL to_c(fcls , cls , 4_C_INT64_T)
CALL to_c(ftext, text, 45_C_INT64_T)
CALL to_c(fmsgs, msgs, 560_C_INT64_T)
!
END SUBROUTINE capi_last_error
!
!*******************************************************************************
!
SUBROUTINE capi_seed(value) BIND(C, NAME='dc_seed')
INTEGER(C_INT), VALUE :: value
CALL dc_seed(INT(value))
END SUBROUTINE capi_seed
!
!*******************************************************************************
!
FUNCTION capi_random_size() BIND(C, NAME='dc_random_size') RESULT(n)
INTEGER(C_INT64_T) :: n
n = dc_random_size()
END FUNCTION capi_random_size
!
!*******************************************************************************
!
FUNCTION capi_get_random_state(state, n) BIND(C, NAME='dc_get_random_state') RESULT(ok)
!
!  Fill state(n) with the generator state; n must be dc_random_size().
!  Returns 0, or -1 for a wrong n.
!
INTEGER(C_INT64_T),                VALUE       :: n
INTEGER(C_INT), DIMENSION(n),      INTENT(OUT) :: state
INTEGER(C_INT) :: ok
INTEGER, DIMENSION(:), ALLOCATABLE :: s
!
ok = -1
state = 0
IF(n /= dc_random_size()) RETURN
ALLOCATE(s(n))
CALL dc_get_random_state(s, INT(n))
state = s
ok = 0
!
END FUNCTION capi_get_random_state
!
!*******************************************************************************
!
FUNCTION capi_set_random_state(state, n) BIND(C, NAME='dc_set_random_state') RESULT(ok)
!
!  Restore a state from dc_get_random_state. Returns 0, or -1 for a wrong n.
!
INTEGER(C_INT64_T),                VALUE      :: n
INTEGER(C_INT), DIMENSION(n),      INTENT(IN) :: state
INTEGER(C_INT) :: ok
INTEGER, DIMENSION(:), ALLOCATABLE :: s
!
ok = -1
IF(n /= dc_random_size()) RETURN
ALLOCATE(s(n))
s = state
CALL dc_set_random_state(s, INT(n))
ok = 0
!
END FUNCTION capi_set_random_state
!
!*******************************************************************************
!
FUNCTION capi_data_sets() BIND(C, NAME='dc_data_sets') RESULT(n)
!  Number of KUPLOT data sets loaded
INTEGER(C_INT) :: n
n = get_data_sets()
END FUNCTION capi_data_sets
!
!*******************************************************************************
!
FUNCTION capi_data_length(ik) BIND(C, NAME='dc_data_length') RESULT(n)
!
!  Number of points of KUPLOT data set ik (1-based); a grid has nx*ny.
!  -1 if there is no such data set.
!
INTEGER(C_INT), VALUE :: ik
INTEGER(C_INT64_T) :: n
!
n = -1
IF(ik < 1 .OR. ik > get_data_sets()) RETURN
n = get_data_length(INT(ik))
!
END FUNCTION capi_data_length
!
!*******************************************************************************
!
FUNCTION capi_get_data(ik, x, y, n) BIND(C, NAME='dc_get_data') RESULT(ok)
!
!  x and y of KUPLOT data set ik; n must be dc_data_length(ik).
!  Returns 0, or -1 (no such data set, or a wrong n).
!
INTEGER(C_INT),                    VALUE       :: ik
INTEGER(C_INT64_T),                VALUE       :: n
REAL(C_FLOAT), DIMENSION(n),       INTENT(OUT) :: x
REAL(C_FLOAT), DIMENSION(n),       INTENT(OUT) :: y
INTEGER(C_INT) :: ok
REAL, DIMENSION(:), ALLOCATABLE :: fx, fy
!
ok = -1
x = 0.0
y = 0.0
IF(capi_data_length(ik) /= n) RETURN
ALLOCATE(fx(n), fy(n))
CALL get_data(INT(ik), fx, fy, INT(n))
x = fx
y = fy
ok = 0
!
END FUNCTION capi_get_data
!
!*******************************************************************************
!
FUNCTION capi_get_data_2d(ik, x, y, z, n, nx, ny) BIND(C, NAME='dc_get_data_2d') RESULT(ok)
!
!  A KUPLOT grid data set: x(n), y(n), z(n) as get_data_3d in
!  diffuse_suite.f90 (x and y per point, z(ix, iy) with y running fastest), and the
!  grid size nx, ny. n must be dc_data_length(ik). Returns 0, or -1 (no
!  such data set, a wrong n, or not a grid).
!
USE kuplot_mod, ONLY: lni           ! only this: KUPLOT has its own x, y, z, nx, ny
!
INTEGER(C_INT),                    VALUE       :: ik
INTEGER(C_INT64_T),                VALUE       :: n
REAL(C_FLOAT), DIMENSION(n),       INTENT(OUT) :: x, y, z
INTEGER(C_INT64_T),                INTENT(OUT) :: nx, ny
INTEGER(C_INT) :: ok
REAL, DIMENSION(:), ALLOCATABLE :: fx, fy, fz
INTEGER, DIMENSION(1) :: npx, npy
!
ok = -1
x = 0.0
y = 0.0
z = 0.0
nx = 0
ny = 0
IF(capi_data_length(ik) /= n) RETURN
IF(.NOT. lni(ik)) RETURN                 ! not a grid
ALLOCATE(fx(n), fy(n), fz(n))
CALL get_data_3d(INT(ik), fx, fy, fz, npx, npy, INT(n))
x = fx
y = fy
z = fz
nx = npx(1)
ny = npy(1)
ok = 0
!
END FUNCTION capi_get_data_2d
!
!*******************************************************************************
!  The crystal DISCUS holds (crystal_mod), as arrays (WP3, first slice: out)
!*******************************************************************************
!
SUBROUTINE capi_crystal_size(natoms, ntypes) BIND(C, NAME='dc_crystal_size')
!
!  Number of atoms in the crystal and of atom types (types run 0..ntypes;
!  type 0 is a vacancy, "VOID").
!
USE crystal_mod, ONLY: cr_natoms, cr_nscat
INTEGER(C_INT64_T), INTENT(OUT) :: natoms
INTEGER(C_INT),     INTENT(OUT) :: ntypes
!
natoms = cr_natoms
ntypes = cr_nscat
!
END SUBROUTINE capi_crystal_size
!
!*******************************************************************************
!
SUBROUTINE capi_get_cell(lattice, spcgr, icc, ncatoms) BIND(C, NAME='dc_get_cell')
!
!  Lattice (a, b, c in A, alpha, beta, gamma in degrees), space group symbol
!  (16 characters, blank padded), unit cells along a, b, c and atoms per
!  unit cell.
!
USE crystal_mod, ONLY: cr_a0, cr_win, cr_spcgr, cr_icc, cr_ncatoms
REAL(C_DOUBLE),         DIMENSION(6),  INTENT(OUT) :: lattice
CHARACTER(KIND=C_CHAR), DIMENSION(16), INTENT(OUT) :: spcgr
INTEGER(C_INT64_T),     DIMENSION(3),  INTENT(OUT) :: icc
INTEGER(C_INT64_T),                    INTENT(OUT) :: ncatoms
!
lattice(1:3) = cr_a0
lattice(4:6) = cr_win
CALL to_c(cr_spcgr, spcgr, 16_C_INT64_T)
icc = cr_icc
ncatoms = cr_ncatoms
!
END SUBROUTINE capi_get_cell
!
!*******************************************************************************
!
FUNCTION capi_get_atoms(n, pos, types, prop, mole) BIND(C, NAME='dc_get_atoms') RESULT(ok)
!
!  The atoms: fractional coordinates pos(3, n) (in C/NumPy order: n rows of
!  x, y, z), atom type (0..ntypes), property flags and molecule number
!  (0: none). n must be the number of atoms. Returns 0, or -1 for a wrong n.
!
USE crystal_mod, ONLY: cr_natoms, cr_pos, cr_iscat, cr_prop, cr_mole
INTEGER(C_INT64_T),                     VALUE       :: n
REAL(C_DOUBLE),     DIMENSION(3, n),    INTENT(OUT) :: pos
INTEGER(C_INT),     DIMENSION(n),       INTENT(OUT) :: types, prop, mole
INTEGER(C_INT) :: ok
INTEGER(C_INT64_T) :: i
!
ok = -1
pos = 0.0_C_DOUBLE
types = 0
prop = 0
mole = 0
IF(n /= cr_natoms) RETURN
ok = 0
IF(n == 0) RETURN
IF(ALLOCATED(cr_pos))   pos = cr_pos(:, 1:n)
IF(ALLOCATED(cr_iscat)) THEN
   DO i = 1, n
      types(i) = cr_iscat(1, i)
   ENDDO
ENDIF
IF(ALLOCATED(cr_prop))  prop = cr_prop(1:n)
IF(ALLOCATED(cr_mole))  mole = cr_mole(1:n)
!
END FUNCTION capi_get_atoms
!
!*******************************************************************************
!
FUNCTION capi_get_atom_types(ntypes, names, biso, occ) BIND(C, NAME='dc_get_atom_types') RESULT(ok)
!
!  The atom types 0..ntypes: names (4 characters each, blank padded, one
!  block of 4*(ntypes+1) characters), isotropic displacement parameter B
!  (A^2) and occupancy. ntypes must be the number from dc_crystal_size.
!  Returns 0, or -1 for a wrong ntypes.
!
USE crystal_mod, ONLY: cr_nscat, cr_at_lis, cr_dw, cr_occ
INTEGER(C_INT),                                 VALUE       :: ntypes
CHARACTER(KIND=C_CHAR), DIMENSION(4*(ntypes+1)), INTENT(OUT) :: names
REAL(C_DOUBLE),         DIMENSION(0:ntypes),    INTENT(OUT) :: biso, occ
INTEGER(C_INT) :: ok
INTEGER :: i, k
!
ok = -1
names = ' '
biso = 0.0_C_DOUBLE
occ = 0.0_C_DOUBLE
IF(ntypes /= cr_nscat) RETURN
IF(.NOT. (ALLOCATED(cr_at_lis) .AND. ALLOCATED(cr_dw) .AND. ALLOCATED(cr_occ))) RETURN
DO i = 0, ntypes
   DO k = 1, 4
      names(4*i + k) = cr_at_lis(i)(k:k)
   ENDDO
   biso(i) = cr_dw(i)
   occ(i)  = cr_occ(i)
ENDDO
ok = 0
!
END FUNCTION capi_get_atom_types
!
!*******************************************************************************
!  A crystal from arrays (WP3, second slice: in)
!*******************************************************************************
!
FUNCTION capi_set_crystal(lattice, spcgr, nsp, icc, ncatoms, ntypes, names, biso, occ, &
                          has_anis, anis, natoms, pos, types, prop) BIND(C, NAME='dc_set_crystal') RESULT(ok)
!
!  Replace the crystal DISCUS holds, as "read stru" does for a file, from
!  arrays: lattice (a, b, c, alpha, beta, gamma), space group symbol
!  (recorded, not applied: pos holds every atom), unit cells along a, b, c
!  and atoms per unit cell, the atom types 1..ntypes (names: 4 characters
!  each in one block, B in A^2, occupancy; if has_anis /= 0 also the
!  anisotropic ADPs anis(6, ntypes): U11, U22, U33, U23, U13, U12 in A^2, as
!  in a CIF; one value U11 > 0 with the rest 0 is an isotropic U) and the
!  atoms (pos(3, natoms)
!  fractional, i.e. C/NumPy order n rows of x, y, z; type 1..ntypes;
!  property flags).
!
!  The steps follow do_readstru / do_readstru_disk / readstru and the read
!  menu in discus/prog/structur.f90: reset, allocate, header, atoms, cell
!  counts, then lattice, symmetry matrices and ADPs. Returns 0, or -1 (an
!  error is recorded for dc_last_error: unknown space group, a type out of
!  range ...).
!
USE crystal_mod
USE chem_mod,       ONLY: chem_period, chem_quick, chem_purge
USE diffuse_mod,    ONLY: four_last, FOUR_NN
USE discus_allocate_appl_mod, ONLY: alloc_crystal_scat, alloc_crystal_nmax, alloc_unitcell
USE discus_config_mod, ONLY: NMAX, MAXSCAT
USE prep_anis_mod,  ONLY: prep_anis
USE spcgr_apply,    ONLY: setup_lattice, get_symmetry_matrices
USE stack_rese_mod, ONLY: do_stack_rese
USE structur,       ONLY: rese_cr, spcgr_no
USE errlist_mod,    ONLY: ier_num, ier_typ, ier_msg, ER_APPL
USE lib_errlist_func, ONLY: errlist
USE precision_mod,  ONLY: PREC_DP
!
REAL(C_DOUBLE),         DIMENSION(6),            INTENT(IN) :: lattice
INTEGER(C_INT64_T),                     VALUE                :: nsp
CHARACTER(KIND=C_CHAR), DIMENSION(nsp),          INTENT(IN) :: spcgr
INTEGER(C_INT64_T),     DIMENSION(3),            INTENT(IN) :: icc
INTEGER(C_INT64_T),                     VALUE                :: ncatoms
INTEGER(C_INT),                         VALUE                :: ntypes
CHARACTER(KIND=C_CHAR), DIMENSION(4*ntypes),     INTENT(IN) :: names
REAL(C_DOUBLE),         DIMENSION(ntypes),       INTENT(IN) :: biso, occ
INTEGER(C_INT),                         VALUE                :: has_anis
REAL(C_DOUBLE),         DIMENSION(6, ntypes),    INTENT(IN) :: anis
INTEGER(C_INT64_T),                     VALUE                :: natoms
REAL(C_DOUBLE),         DIMENSION(3, natoms),    INTENT(IN) :: pos
INTEGER(C_INT),         DIMENSION(natoms),       INTENT(IN) :: types, prop
INTEGER(C_INT) :: ok
!
INTEGER, PARAMETER :: MAXW = 1
REAL(KIND=PREC_DP), DIMENSION(MAXW) :: werte
INTEGER :: ianz, i, k, j, ncheck, ncatoms_saved
INTEGER, DIMENSION(:), ALLOCATABLE :: type_adp
LOGICAL :: lout
!
ok = -1
CALL set_it
IF(ier_num /= 0) THEN
   CALL errlist                      ! records the error for dc_last_error
   RETURN
ENDIF
ok = 0
!
CONTAINS
!
SUBROUTINE set_it
IF(ntypes < 1 .OR. natoms < 1 .OR. ncatoms < 1 .OR. MINVAL(icc) < 1) THEN
   ier_num = -6
   ier_typ = ER_APPL
   ier_msg(1) = 'dc_set_crystal: empty crystal or bad cell counts'
   RETURN
ENDIF
IF(MINVAL(types) < 1 .OR. MAXVAL(types) > ntypes) THEN
   ier_num = -6
   ier_typ = ER_APPL
   ier_msg(1) = 'dc_set_crystal: atom type out of range'
   RETURN
ENDIF
!
CALL rese_cr
IF(natoms > NMAX .OR. ntypes > MAXSCAT) THEN
   CALL alloc_crystal_scat(MAX(INT(ntypes) + 2, MAXSCAT))
   CALL alloc_crystal_nmax(MAX(INT(natoms) + 10, NMAX))
   IF(ier_num /= 0) RETURN
ENDIF
!
!  Header: name, space group, lattice, atom types
cr_name  = 'from Python'
cr_spcgr = from_c(spcgr, MIN(nsp, 16_C_INT64_T))
cr_set   = 'abc'
cr_a0    = lattice(1:3)
cr_win   = lattice(4:6)
cr_nscat = ntypes
cr_at_lis(0) = 'VOID'
cr_dw(0)  = 0.0_PREC_DP
cr_occ(0) = 1.0_PREC_DP
DO i = 1, ntypes
   cr_at_lis(i) = ' '
   DO k = 1, 4
      cr_at_lis(i)(k:k) = names(4*(i-1) + k)
   ENDDO
   cr_dw(i)  = biso(i)
   cr_occ(i) = occ(i)
ENDDO
cr_anis(:, 0:ntypes) = 0.0_PREC_DP
IF(has_anis /= 0) cr_anis(:, 1:ntypes) = anis
werte = 0.0_PREC_DP
ianz = 0
CALL spcgr_no(ianz, MAXW, werte)               ! space group number, crystal system
IF(ier_num /= 0) THEN
   ier_msg(1) = 'Space group: ' // cr_spcgr
   RETURN
ENDIF
!
!  Atoms, as struc_read_atoms sets them
cr_natoms = INT(natoms)
cr_dim(:, 1) =  1.0D10
cr_dim(:, 2) = -1.0D10
as_natoms = 0
DO i = 1, cr_natoms
   cr_pos(:, i)   = pos(:, i)
   cr_iscat(1, i) = types(i)
   cr_iscat(2:, i) = 1
   cr_ianis(i)    = types(i)
   cr_prop(i)     = prop(i)
   cr_mole(i)     = 0
   cr_surf(:, i)  = 0
   cr_magn(:, i)  = 0.0_PREC_DP
   cr_valu(i)     = 0.0_PREC_DP
   DO j = 1, 3
      cr_dim(j, 1) = MIN(cr_dim(j, 1), cr_pos(j, i))
      cr_dim(j, 2) = MAX(cr_dim(j, 2), cr_pos(j, i))
   ENDDO
ENDDO
!  First atom of each type inside the first unit cell: the "asymmetric unit"
!  list that struc_read_atoms keeps
DO k = 1, ntypes
   DO i = 1, cr_natoms
      IF(types(i) == k .AND. ALL(pos(:, i) >= 0.0_C_DOUBLE) .AND. ALL(pos(:, i) < 1.0_C_DOUBLE)) THEN
         as_natoms = as_natoms + 1
         as_at_lis(as_natoms) = cr_at_lis(k)
         as_pos(:, as_natoms) = cr_pos(:, i)
         as_iscat(as_natoms)  = k
         as_prop(as_natoms)   = cr_prop(i)
         as_dw(as_natoms)     = cr_dw(k)
         EXIT
      ENDIF
   ENDDO
ENDDO
DO j = 1, 3
   cr_dim0(j, 1) = REAL(NINT(cr_dim(j, 1)))
   cr_dim0(j, 2) = REAL(NINT(cr_dim(j, 2)))
ENDDO
!
!  Crystal size, as do_readstru_disk with "ncell" in the file
cr_icc     = INT(icc)
cr_ncatoms = INT(ncatoms)
cr_ncreal  = INT(ncatoms)
IF(cr_natoms /= cr_icc(1)*cr_icc(2)*cr_icc(3)*cr_ncatoms) THEN
   chem_period(:) = .FALSE.
   chem_quick     = .FALSE.
ELSE
   chem_period(:) = .TRUE.
   chem_quick     = .TRUE.
ENDIF
CALL alloc_unitcell(cr_ncatoms)
cr_is_sym = 1
chem_purge = .FALSE.
!
!  As the read menu after a structure: lattice, symmetry, ADPs
lout = .FALSE.
CALL setup_lattice(cr_a0, cr_ar, cr_eps, cr_gten, cr_reps, cr_rten, cr_win, cr_wrez, &
                   cr_v, cr_vr, lout, cr_gmat, cr_fmat, cr_cartesian,                 &
                   cr_tran_g, cr_tran_gi, cr_tran_f, cr_tran_fi)
CALL get_symmetry_matrices
IF(ier_num /= 0) RETURN
!  ADPs. prep_anis sets the tensor of the first 5 * cr_ncatoms atoms only.
!  Here every atom is "first" (cr_iscat(2,i) = 1, nothing is rotated), so
!  all atoms of a type share one tensor: copy it to the others. If a type
!  first appears beyond that range, let prep_anis look at every atom.
ncheck = MIN(cr_natoms, 5 * cr_ncatoms)
ALLOCATE(type_adp(0:cr_nscat))
type_adp = 0
DO i = 1, ncheck
   type_adp(cr_iscat(1, i)) = 1
ENDDO
ncatoms_saved = cr_ncatoms
IF(ANY(type_adp(1:cr_nscat) == 0)) cr_ncatoms = cr_natoms
CALL prep_anis(cr_natoms, .TRUE.)
cr_ncatoms = ncatoms_saved
IF(ier_num /= 0) RETURN
type_adp = 0
DO i = 1, cr_natoms
   IF(type_adp(cr_iscat(1, i)) == 0) type_adp(cr_iscat(1, i)) = cr_iscat(3, i)
   IF(i > ncheck) cr_iscat(3, i) = type_adp(cr_iscat(1, i))
ENDDO
DEALLOCATE(type_adp)
CALL do_stack_rese
four_last = FOUR_NN
END SUBROUTINE set_it
!
END FUNCTION capi_set_crystal
!
!*******************************************************************************
!  Displacement parameters (WP3, third slice)
!*******************************************************************************
!
FUNCTION capi_get_type_anis(ntypes, anis) BIND(C, NAME='dc_get_type_anis') RESULT(ok)
!
!  The ADPs of the atom types 0..ntypes as given (file or dc_set_crystal):
!  anis(6, 0:ntypes) = U11, U22, U33, U23, U13, U12 (A^2; all 0: isotropic,
!  from B). Returns 0, or -1 for a wrong ntypes.
!
USE crystal_mod, ONLY: cr_nscat, cr_anis
INTEGER(C_INT),                         VALUE       :: ntypes
REAL(C_DOUBLE), DIMENSION(6, 0:ntypes), INTENT(OUT) :: anis
INTEGER(C_INT) :: ok
!
ok = -1
anis = 0.0_C_DOUBLE
IF(ntypes /= cr_nscat .OR. .NOT. ALLOCATED(cr_anis)) RETURN
anis = cr_anis(:, 0:ntypes)
ok = 0
!
END FUNCTION capi_get_type_anis
!
!*******************************************************************************
!
FUNCTION capi_adp_count() BIND(C, NAME='dc_adp_count') RESULT(n)
!
!  Number of distinct displacement tensors DISCUS made from the types and the
!  symmetry (prep_anis): the list dc_get_adps returns.
!
USE crystal_mod, ONLY: cr_nanis
INTEGER(C_INT64_T) :: n
n = cr_nanis
END FUNCTION capi_adp_count
!
!*******************************************************************************
!
FUNCTION capi_get_adps(nanis, uij, prin) BIND(C, NAME='dc_get_adps') RESULT(ok)
!
!  The distinct displacement tensors: uij(6, nanis) = U11, U22, U33, U23,
!  U13, U12 (A^2, CIF convention) and the principal axes prin(4, 3, nanis):
!  for axis k, prin(1:3, k, j) the unit vector (Cartesian) and prin(4, k, j)
!  the mean-square displacement along it. Returns 0, or -1 for a wrong nanis.
!
USE crystal_mod, ONLY: cr_nanis, cr_anis_full, cr_prin
INTEGER(C_INT64_T),                         VALUE       :: nanis
REAL(C_DOUBLE),     DIMENSION(6, nanis),    INTENT(OUT) :: uij
REAL(C_DOUBLE),     DIMENSION(4, 3, nanis), INTENT(OUT) :: prin
INTEGER(C_INT) :: ok
!
ok = -1
uij = 0.0_C_DOUBLE
prin = 0.0_C_DOUBLE
IF(nanis /= cr_nanis) RETURN
ok = 0
IF(nanis == 0) RETURN
IF(ALLOCATED(cr_anis_full)) uij  = cr_anis_full(:, 1:nanis)
IF(ALLOCATED(cr_prin))      prin = cr_prin(:, :, 1:nanis)
!
END FUNCTION capi_get_adps
!
!*******************************************************************************
!
FUNCTION capi_get_atom_adps(n, adp) BIND(C, NAME='dc_get_atom_adps') RESULT(ok)
!
!  For each atom, its displacement tensor: index 1..nanis into dc_get_adps.
!  n must be the number of atoms. Returns 0, or -1 for a wrong n.
!
USE crystal_mod, ONLY: cr_natoms, cr_iscat
INTEGER(C_INT64_T),               VALUE       :: n
INTEGER(C_INT),     DIMENSION(n), INTENT(OUT) :: adp
INTEGER(C_INT) :: ok
INTEGER(C_INT64_T) :: i
!
ok = -1
adp = 0
IF(n /= cr_natoms) RETURN
ok = 0
IF(n == 0 .OR. .NOT. ALLOCATED(cr_iscat)) RETURN
DO i = 1, n
   adp(i) = cr_iscat(3, i)
ENDDO
!
END FUNCTION capi_get_atom_adps
!
!*******************************************************************************
!  Molecules (WP3, third slice)
!*******************************************************************************
!
SUBROUTINE capi_molecule_size(nmole, ntypes, natoms) BIND(C, NAME='dc_molecule_size')
!
!  Number of molecules, of molecule types and of atoms in all molecules.
!
USE molecule_mod, ONLY: mole_num_mole, mole_num_type, mole_off, mole_len
INTEGER(C_INT64_T), INTENT(OUT) :: nmole
INTEGER(C_INT),     INTENT(OUT) :: ntypes
INTEGER(C_INT64_T), INTENT(OUT) :: natoms
!
nmole  = mole_num_mole
ntypes = mole_num_type
natoms = 0
IF(mole_num_mole > 0) natoms = mole_off(mole_num_mole) + mole_len(mole_num_mole)
!
END SUBROUTINE capi_molecule_size
!
!*******************************************************************************
!
FUNCTION capi_get_molecules(nmole, natoms, mtype, mchar, mlen, atoms) &
         BIND(C, NAME='dc_get_molecules') RESULT(ok)
!
!  The molecules: type, character (0: atoms; 1 cube, 2 cylinder, 3 sphere,
!  4 edge; negative: domains), number of atoms, and the atoms of all
!  molecules one after the other (1-based atom numbers; molecule m has
!  mlen(m) atoms). Returns 0, or -1 for wrong sizes.
!
USE molecule_mod, ONLY: mole_num_mole, mole_type, mole_char, mole_off, mole_len, mole_cont
INTEGER(C_INT64_T),                    VALUE       :: nmole, natoms
INTEGER(C_INT),     DIMENSION(nmole),  INTENT(OUT) :: mtype, mchar, mlen
INTEGER(C_INT),     DIMENSION(natoms), INTENT(OUT) :: atoms
INTEGER(C_INT) :: ok
INTEGER :: m, k, i
!
ok = -1
mtype = 0
mchar = 0
mlen = 0
atoms = 0
IF(nmole /= mole_num_mole) RETURN
IF(nmole > 0) THEN
   IF(natoms /= mole_off(mole_num_mole) + mole_len(mole_num_mole)) RETURN
ELSEIF(natoms /= 0) THEN
   RETURN
ENDIF
k = 0
DO m = 1, mole_num_mole
   mtype(m) = mole_type(m)
   mchar(m) = mole_char(m)
   mlen(m)  = mole_len(m)
   DO i = 1, mole_len(m)
      k = k + 1
      atoms(k) = mole_cont(mole_off(m) + i)
   ENDDO
ENDDO
ok = 0
!
END FUNCTION capi_get_molecules
!
!*******************************************************************************
!
FUNCTION capi_set_molecules(nmole, natoms, mtype, mlen, atoms) &
         BIND(C, NAME='dc_set_molecules') RESULT(ok)
!
!  Group atoms of the crystal DISCUS holds into molecules (after
!  dc_set_crystal): molecule m has type mtype(m) (1, 2, ...) and the
!  mlen(m) atoms that follow in atoms (1-based atom numbers). Replaces all
!  molecules; character "atoms" (0). Returns 0, or -1 (error recorded).
!
USE crystal_mod,   ONLY: cr_natoms, cr_prop, cr_mole
USE molecule_mod
USE prop_para_mod, ONLY: PROP_MOLECULE
USE discus_allocate_appl_mod, ONLY: alloc_molecule
USE errlist_mod,   ONLY: ier_num, ier_typ, ier_msg, ER_APPL
USE lib_errlist_func, ONLY: errlist
USE precision_mod, ONLY: PREC_DP
INTEGER(C_INT64_T),                    VALUE      :: nmole, natoms
INTEGER(C_INT),     DIMENSION(nmole),  INTENT(IN) :: mtype, mlen
INTEGER(C_INT),     DIMENSION(natoms), INTENT(IN) :: atoms
INTEGER(C_INT) :: ok
INTEGER :: m, i, k, ntype
!
ok = -1
IF(SUM(INT(mlen, C_INT64_T)) /= natoms .OR. ANY(mlen < 1) .OR. ANY(mtype < 1) .OR. &
   ANY(atoms < 1) .OR. ANY(atoms > cr_natoms)) THEN
   ier_num = -6
   ier_typ = ER_APPL
   ier_msg(1) = 'dc_set_molecules: sizes, types or atom numbers out of range'
   CALL errlist
   RETURN
ENDIF
!  Forget the old molecules
DO i = 1, cr_natoms
   cr_mole(i) = 0
   cr_prop(i) = IBCLR(cr_prop(i), PROP_MOLECULE)
ENDDO
ntype = 0
IF(nmole > 0) ntype = MAXVAL(mtype)
IF(nmole > MOLE_MAX_MOLE .OR. ntype > MOLE_MAX_TYPE .OR. natoms > MOLE_MAX_ATOM) THEN
   CALL alloc_molecule(MAX(1, MOLE_MAX_GENE), MAX(1, MOLE_MAX_SYMM), MAX(INT(nmole) + 1, MOLE_MAX_MOLE), &
                       MAX(ntype + 1, MOLE_MAX_TYPE), MAX(INT(natoms) + 1, MOLE_MAX_ATOM))
   IF(ier_num /= 0) THEN
      CALL errlist
      RETURN
   ENDIF
ENDIF
mole_num_mole = INT(nmole)
mole_num_type = ntype
mole_num_curr = 0
k = 0
DO m = 1, mole_num_mole
   mole_type(m) = mtype(m)
   mole_char(m) = MOLE_ATOM
   mole_dens(m) = 0.0_PREC_DP
   mole_file(m) = ' '
   mole_fuzzy(m) = 0.0_PREC_DP
   mole_off(m)  = k
   mole_len(m)  = mlen(m)
   DO i = 1, mlen(m)
      k = k + 1
      mole_cont(k) = atoms(k)
      cr_mole(atoms(k)) = m
      cr_prop(atoms(k)) = IBSET(cr_prop(atoms(k)), PROP_MOLECULE)
   ENDDO
ENDDO
mole_num_atom = k
IF(ntype > 0) THEN
   mole_biso(1:ntype) = 0.0_PREC_DP
   mole_clin(1:ntype) = 0.0_PREC_DP
   mole_cqua(1:ntype) = 0.0_PREC_DP
ENDIF
ok = 0
!
END FUNCTION capi_set_molecules
!
END MODULE diffuse_capi
