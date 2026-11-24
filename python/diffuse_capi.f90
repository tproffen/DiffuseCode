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
!  (the module "suite" in suite_python.f90 already has dc_init etc.).
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
INTEGER(C_INT), PARAMETER :: API_VERSION = 1
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
!  A KUPLOT grid data set: x(n), y(n), z(n) as in get_data_3d of the f2py
!  module (x and y per point, z(ix, iy) with y running fastest), and the
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
END MODULE diffuse_capi
