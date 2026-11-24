# WP2 – Binding decision: f2py or a C interface

WP2 step 5 (migration plan, section 7; survey in `WP2_Core_Survey.md`).
Written Nov 2026, after steps 1–4; option E (replace the Fortran) added
later. **Status: proposal; the decision is the owner's.**

## The question

Python reaches DISCUS through `suite_python`, an extension module that
f2py generates from `python/suite_python.f90`. The plan foresaw f2py for
the start and a `bind(C)` interface later. The question is whether to
switch, to what, and when – before WP3 (crystal data model) adds many new
functions with arrays in both directions.

A further option was raised later: replace the Fortran itself (C++, or
Python with NumPy/Numba). It is discussed in
[Option E](#option-e-replace-the-fortran-c-or-python-with-numpy-numba),
with an inventory of how much Fortran would actually have to be replicated.

## What the Python side needs from the binding

From the work so far and from WP3/WP4:

* coarse calls: run commands, read a value, read and write arrays
  (crystal: positions, types, occupancies, ADPs; results: maps,
  patterns), error records, random state;
* arrays of a size known only on the Fortran side (two-step "size, then
  fill" works);
* the Python GIL released during long calls (live output, responsive
  notebooks);
* simple builds for users without admin rights and, later, wheels for
  Linux and macOS (WP1).

Per-call overhead does not matter: a call runs a DISCUS command or moves a
whole array.

## What f2py has cost so far

Concrete items from this migration, not general criticism:

| Issue | Where it showed | Workaround |
|---|---|---|
| f2py crashed (`StopIteration`) on Ubuntu 22.04 with a user-installed NumPy | the owner's WSL2 build | `cmake/FindF2PY.cmake` runs `python -m numpy.f2py` (PR #39) |
| Python's `distutils` missing in pip's isolated build | `pip install .` (WP1) | `EXT_SUFFIX` from `sysconfig` instead of FindPythonExtensions |
| Hundreds of lines of f2py messages per build | every build | `--quiet` |
| `intent(out)` ignored for module procedures | `dc_last_error`, `dc_get_random_state` (WP2) | `!f2py intent(out)` directives |
| Strings come back as padded `bytes`; no arrays of strings | error details (7 × 80 characters joined into one string) | decode and split in Python |
| Releasing the GIL needs a directive | live output (WP2 step 1) | `!f2py threadsafe` |
| The extension depends on the Python version and is built against NumPy's C API | packaging | NumPy at build time; one wheel per Python version (cp310 … cp313) and platform |

Each was solvable, but they cluster in the build and packaging, which is
where users without admin rights get stuck.

## Options

| | A: keep f2py | B: `bind(C)` + ctypes | C: `bind(C)` + nanobind / pybind11 | D: `bind(C)` + Cython |
|---|---|---|---|---|
| Fortran side | as now | small C-API layer (`bind(C)`, C types only) | same as B | same as B |
| Python side | generated | one Python file declaring the functions (`ctypes`, `numpy.ctypeslib`) | C++ glue, compiled | `.pyx` glue, compiled |
| Build needs | NumPy + f2py at build time | none beyond the Fortran build | C++ compiler, nanobind/pybind11 | Cython, C compiler |
| Depends on Python version | yes (extension module) | **no** (plain shared library) | yes (stable ABI only from Python 3.12, nanobind) | yes |
| Wheels | per Python version and platform | **one per platform** (`py3-none-…`) | per version (or 3.12+ stable ABI) | per version |
| GIL during calls | released with `!f2py threadsafe` | **released automatically** | explicit | explicit (`with nogil`) |
| Arrays | NumPy arrays, caller allocates | NumPy arrays via `ndpointer`, caller allocates | zero-copy views, rich types | typed memoryviews |
| Wrong signature | f2py checks | **crash, not an exception** (no check) | compile-time check | compile-time check |
| Other languages (C, Julia, a service layer in WP10) | no | **yes, the same C API** | via the C API | via the C API |
| Effort to switch | none | 2–4 days (below) | 1–2 weeks | about a week |

## Proof of concept (option B)

To test B before recommending it, a small `bind(C)` layer over the existing
core functions (`dc_init`, running commands from memory, `get_value`,
`dc_last_error`, KUPLOT data) was linked against the normal build into
`libdiffuse_core.so` (about 8 MB, the same libraries as the f2py module)
and called from Python with `ctypes` only:

* commands run, values come back (`r[1] = 6*7` → 42.0), the error record
  arrives (−8, "Unknown command"), KUPLOT data fills NumPy arrays;
* another Python thread kept running during a 1-second DISCUS call
  (10 of 10 ticks): ctypes releases the GIL by itself;
* the library has **no undefined Python symbols** (`nm -D`): it does not
  depend on the Python version or on NumPy's C API.

Each C function was 3–10 lines of Fortran: copy a `c_char` array to a
Fortran string, call the existing routine, pass plain numbers back. The
proof of concept is not committed; step 1 below would build the real
layer.

## Recommendation: B, a C interface with ctypes, in steps

B removes the build and packaging problems above, gives one wheel per
platform instead of one per Python version, releases the GIL without
directives, and produces a C interface that other languages and the WP10
service layer can use as well. Its weakness, unchecked signatures, is
contained: all declarations live in one Python file, and every function is
covered by tests (a wrong declaration crashes the test run immediately).
The richer types of C (nanobind) are not needed for coarse calls; they can
be added on top of the same C interface if WP3 ever needs them.

Steps, each with green CI before the next:

1. **C interface.** `python/diffuse_capi.f90`: `bind(C)` functions for
   everything `suite_python` offers now (prefix `dc_`, C types only,
   strings as `char` + length, arrays with explicit size). CMake builds it
   as a shared library `libdiffuse_core` next to the f2py module and
   installs it into the package. About 1 day.
   *Done (Nov 2026):* `python/diffuse_capi.f90`, header
   `python/diffuse_core.h`, tests `tests/python/test_capi.py` (every
   function, called with ctypes in a separate process). Two differences from
   `suite_python`: `dc_value` returns a double (about 12 significant
   digits, as many as the expression evaluator keeps; f2py gave single
   precision), and `dc_set_value`, like `set_value`, sets user variables
   only (`r[]` and other built-in arrays are set with a command).
2. **Python loader.** `diffuse/_core.py` loads the library with ctypes and
   declares every function in one place; `_session.py` uses it when it is
   there and falls back to `suite_python` (as it does for older module
   builds today). The whole test suite runs on both. About 1 day.
   *Done (Nov 2026):* `CtypesCore` and `F2pyCore` in `_core.py` offer the
   same methods; `Session` uses `_core.load()`, which takes the C library
   when it is there. `DIFFUSE_BACKEND=ctypes|f2py` forces one; CI runs
   pytest on f2py in the build job and on ctypes in the installed-package
   and conda jobs. `Session.backend` tells which is in use. With ctypes,
   `value()` returns about 12 significant digits instead of 7.
3. **Drop f2py.** Remove `suite_python`, `cmake/FindF2PY.cmake` and the
   NumPy build dependency; `pyproject.toml` builds a platform wheel that
   works for every Python 3 (`wheel.py-api = "py3"` in scikit-build-core).
   About half a day, plus CI.
   *Done (Nov 2026):* `python/CMakeLists.txt` builds only
   `libdiffuse_core`; the f2py, NumPy, Cython and PythonExtensions CMake
   modules are gone, and `DIFFUSE_PYTHON` no longer looks for Python at all.
   `python/suite_python.f90` is now `python/diffuse_suite.f90` (the same
   Fortran module `suite`, without the `!f2py` lines) and is compiled into
   the library only. `F2pyCore` and `DIFFUSE_BACKEND` are removed; the
   build needs no NumPy, and the wheel is `py3-none-<platform>`.
4. **From WP3 on**, new functions (crystal arrays in and out) are written
   directly as `bind(C)` functions.
   *Started (Nov 2026), crystal out:* `dc_crystal_size`, `dc_get_cell`,
   `dc_get_atoms` (positions, types, properties, molecules) and
   `dc_get_atom_types` (names, B, occupancy); C interface version 2.
   Pattern "ask, then fill"; positions as `double[n][3]`, which is the
   memory layout of `cr_pos(3, n)`. `Crystal.atoms()` reads them instead
   of saving and parsing a `.stru` file; new `Crystal.atom_types()` and
   `Crystal.space_group`. Next: crystal in (Python arrays into DISCUS).

## Risks and how to handle them

* **Signature mismatch** (ctypes cannot check): one declaration table, a
  test that calls every declared function, and simple types only.
* **Finding the library**: it sits inside the installed package and is
  loaded by its full path; on macOS the install name must use `@rpath`
  (handled with CMake's `INSTALL_RPATH`).
* **Wheels need libgfortran, FFTW, HDF5 bundled**: the same as with f2py;
  `auditwheel` (Linux) and `delocate` (macOS) bundle the dependencies of
  any shared library in the wheel.
* **Fortran strings and arrays**: `c_char` arrays with an explicit length
  and explicit-size arrays only (no assumed-shape arrays, which would need
  C descriptors); see the next section.

## Strings and arrays across the C interface

### The two constraints

**Strings.** A normal Fortran string argument (`character(len=*)`) carries
its length as an extra hidden argument, passed in a compiler-specific way,
so a `bind(C)` function cannot take one. It receives a C byte array plus an
explicit length and copies it into a Fortran string (`to_fortran` in the
proof of concept):

```fortran
subroutine dc_c_run(text, n) bind(C)
   integer(c_int), value :: n
   character(kind=c_char), dimension(n), intent(in) :: text   ! bytes + length
```

Strings going back work the same way: Python provides a buffer of fixed
size, Fortran fills and pads it.

**Arrays.** Fortran has two kinds of array argument:

* *explicit-size* (`real :: x(n)`, `n` passed alongside): just a pointer to
  the first element. C and ctypes handle it directly; NumPy arrays go in
  with `numpy.ctypeslib.ndpointer`.
* *assumed-shape* (`real :: x(:)`), used by most modern DISCUS routines: the
  compiler passes a hidden descriptor (address, bounds, strides). Since
  Fortran 2018 a `bind(C)` routine may take these too, but the caller then
  has to build a `CFI_cdesc_t` (from `ISO_Fortran_binding.h`). ctypes cannot
  include C headers; the layout would have to be copied by hand and can
  differ between compilers and versions. Not worth the risk.

So the C layer uses explicit-size arrays at the boundary only. Inside the
wrapper nothing changes: passing an explicit-size array on to an existing
routine that expects `x(:)` makes the compiler build the descriptor itself.
**Existing DISCUS code needs no change.**

### What this means later, and how to handle it

1. **Arrays whose size only Fortran knows** (atoms after building a
   crystal, map sizes). Python cannot receive an array Fortran allocates.
   Two patterns:
   * *ask, then fill*: one call returns the size, Python allocates, a second
     call fills (`get_data_length` / `get_data` already work this way).
     Simple and safe; the data is copied.
   * *view Fortran's memory*: Fortran returns the address (`c_loc`) and
     size, NumPy wraps it without copying (`np.ctypeslib.as_array`). Fast
     for large maps, but the view is invalid as soon as DISCUS reallocates
     the array (e.g. the next structure is read). If ever needed: internal
     only, with an explicit copy or a version counter checked before use.

   For WP3/WP4 sizes (up to millions of atoms, maps of a few hundred MB)
   copying is fine; the second pattern matters only if memory gets tight.
2. **Several dimensions.** Fortran stores arrays column by column, NumPy
   row by row by default. Often this works out: `cr_pos(3, natoms)` has
   exactly the memory layout of a NumPy array of shape `(natoms, 3)`. For
   maps the wrapper returns a Fortran-ordered array (`order="F"`) or the
   transposed shape. Easy to get wrong without tests, so every array
   function gets a test with a non-symmetric shape.
3. **Non-contiguous NumPy arrays** (slices such as `a[:, ::2]`) cannot go in
   directly; the Python wrapper converts them (`np.asfortranarray` or
   `np.ascontiguousarray`), at the cost of a copy. In one place: the
   wrapper, not user code.
4. **Arrays of strings** (atom names, `cr_at_lis`, 4 characters each):
   passed as one block of `n × length` bytes and split in Python, as the
   error details are now. Fixed width and padded; a longer name would be
   cut off, so the wrapper checks and raises an error.
5. **Size limits.** `c_int` is 32-bit; counts above about 2.1 billion
   (possible for large 3D grids) need 64-bit (`c_int64_t`) arguments. Use
   64-bit for all sizes from the start.
6. **Derived types.** Fortran types with allocatable arrays, pointers or
   type-bound procedures cannot cross the boundary (molecule lists, the
   macro storage, the `class_*` objects). The C interface offers functions
   per piece of data (get/set positions, types, occupancies, ADPs), not
   structures; the Python `Crystal` class assembles them. More functions,
   each simple.
7. **Logicals and optional arguments.** Fortran `logical` is not C's
   `bool`, and optional arguments across `bind(C)` depend on the compiler.
   Integer flags (0/1) and fixed argument lists.
8. **Callbacks from Fortran into Python** (progress reporting, an "abort"
   check during MMC): possible with `c_funptr` in Fortran and
   `ctypes.CFUNCTYPE` in Python, which takes the GIL back for the callback
   by itself. A Python exception cannot travel through Fortran, so a
   callback signals an abort with its return value and Fortran stops in an
   orderly way. Not needed now (live output works without callbacks); only
   for interactive abort or progress bars.
9. **No signature check.** ctypes trusts the declared types; a mismatch
   crashes instead of raising an error. Safeguards: one declaration table,
   a test per function, and only a few patterns (bytes + length, pointer +
   size, integer flags), so mistakes are easy to spot.

### When this would stop being enough

If a later work package needed rich two-way objects (Python objects living
in Fortran, many small calls, zero-copy everywhere), a C++ layer with
nanobind could sit **on top of the same C interface** and build the
descriptors properly; nothing from option B would be thrown away. What
WP3–WP9 need (whole arrays in and out, coarse commands) does not call for
it, but the route exists.

## Option E: replace the Fortran (C++, or Python with NumPy/Numba)

Raised Nov 2026: instead of binding to the Fortran, recreate the core
calculations in C or C++ (or Python). The question is how much Fortran
would really have to be replicated, since much of it is replaced by the
Python migration anyway.

### Inventory

`tools/fortran_inventory.py` sorts every Fortran file (code lines, without
comments and blank lines) into groups; DISCUS files that hold both a menu
and the physics behind it are split routine by routine (a routine counts as
menu if its name says so or if more than 4 % of its lines call the command
parser). An estimate, good to perhaps ±10 % for the physics core.

| Group | Lines | Share | Future |
|---|---:|---:|---|
| A command language, menus, program frame, macros, expressions | 67 800 | 33 % | replaced by the Python API |
| B plotting (KUPLOT graphics, structure plots) | 10 700 | 5 % | matplotlib, `Crystal.view()` |
| C fitting and optimisation (DIFFEV, REFINE, KUPLOT fit) | 19 500 | 9 % | SciPy / lmfit, Python orchestration |
| D file formats and export | 24 700 | 12 % | Python (h5py, gemmi, NumPy); DISCUS's own formats rewritten in Python |
| E general numerics (FFT, splines, matrices, random, helpers) | 21 300 | 10 % | NumPy / SciPy / FFTW |
| F system (OS, MPI, sockets, timers) | 3 800 | 2 % | Python standard library |
| **G physics core** | **55 500** | **27 %** | **to replicate (or keep in Fortran)** |
| H data tables (scattering factors, profile functions) | 3 700 | 2 % | keep as data (or from xraydb, periodictable) |
| total | 206 900 | | |

About **three quarters of the Fortran goes away with the Python migration
whatever happens** (A–F). The rewrite question concerns the physics core,
about 55 000 lines (45 000–65 000 given the uncertainty), plus the data
tables:

| Physics area | Lines | Notes |
|---|---:|---|
| Fourier, powder, PDF, Patterson, 3D-PDF | 16 800 | well covered by golden tests (WP0, WP4); GPU candidates |
| structure: reading, symmetry, space groups, generators | 14 200 | partly available in libraries (gemmi, spglib: space groups, symmetry operations, settings) |
| molecules, domains, stacking faults, surfaces, waves, modify | 10 400 | specialised, DISCUS-specific |
| Monte Carlo (MMC, RMC) | 8 500 | many energy terms; intertwined with the crystal and chemistry code |
| chemistry: correlations, connectivity, properties | 5 500 | partly expressible with NumPy |

### Effort

Porting validated numerical code, with tests, runs at roughly 100–250 lines
of Fortran a day for a developer who knows the physics. For the 55 000
lines that is about 1–2.5 person-years, before performance work (OpenMP or
GPU) and before extending the tests where WP0 found gaps
(`WP0_Test_Coverage_Gaps.md`). Using libraries for symmetry (gemmi,
spglib) would save part of the structure area. AI-assisted translation
speeds up writing, not validation, which is the larger part.

### Weighing it

For a rewrite: more people know C++ (or Python) than Fortran; GPUs are
more natural; better tools; no Fortran runtime in wheels; a C++ core would
bind to Python with nanobind directly.

Against a rewrite all at once: the physics is decades of validated work,
and every routine would need re-validating; Fortran is good at exactly
this kind of array code; the problems met so far (STOP, global state,
keyboard reads, macro files) were architecture, not language, and are being
fixed – the one serious memory bug found was in C++ (the bundled finufft);
the people who know the physics write Fortran; and a parallel rewrite
delays WP3–WP10 or doubles the work.

### Recommendation: gradual, behind the C interface

The C interface of option B is the seam a gradual move needs. Every core
calculation sits behind a C function; any of them can later be
reimplemented (C++, or Python with NumPy/Numba/JAX), checked against the
golden tests and swapped in without users noticing. f2py offers no such
seam.

1. Option B as proposed (steps 1–3 above).
2. Plan rule: a core area moves out of Fortran only with a concrete reason
   (speed, GPU, maintainability, a missing expert) and only once golden
   tests cover it; otherwise it stays Fortran behind the C interface.
3. A pilot: one well-tested kernel – powder (Debye) or PDF – reimplemented
   both in C++ and in Python with NumPy/Numba, compared with the Fortran
   for results (golden tests), speed and effort. About 1–2 weeks; it gives
   measured numbers instead of the estimates above before anything larger
   is decided. The pilot includes **benchmarks** (owner's request), not
   only equal results: the same calculations timed for the Fortran and
   each rebuild, on one thread and on all cores, for several crystal
   sizes (from a few hundred to millions of atoms) and grid sizes, plus
   peak memory; each run repeated, reporting the median; the machine, the
   compilers and their flags recorded; and kept as a script, so that the
   comparison can be repeated (later also on a GPU).
4. Likely order if the pilot is convincing: Fourier/powder/PDF first
   (best tested, most to gain on GPUs), symmetry next with libraries
   (gemmi/spglib), Monte Carlo and the specialised structure tools last or
   never.

A complete rewrite at once makes sense only with a dedicated team and the
Fortran authors on board – a question of people and funding more than of
code.

## What does not change

The public API (`diffuse.read_structure`, `powder`, `fourier`, `pdf`,
`Crystal.view`, `MacroError`, `seed`, `reset` …), the way commands reach
DISCUS (macro held in memory), output capture and live echo, and all tests.
Only `_session.py`'s lowest layer changes.
