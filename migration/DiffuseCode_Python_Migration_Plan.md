# DiffuseCode: moving to a Python framework

A high-level look at the current Fortran code and options for separating compute from the command language, graphics and I/O, with a web front end as a later goal.

> **Related documents:** progress on WP0 is tracked in [`WP0_Test_Baseline.md`](WP0_Test_Baseline.md); test gaps are in [`WP0_Test_Coverage_Gaps.md`](WP0_Test_Coverage_Gaps.md). `mixscat` is deprecated and not part of the migration.

1.  [What the code is today](#1-what-the-code-is-today)
2.  [Target architecture](#2-target-architecture-shared-by-all-options)
3.  [High-level options](#3-high-level-options)
4.  [Suggested route](#4-suggested-route-b-then-d)
5.  [Main risks and decisions to make early](#5-main-risks-and-decisions-to-make-early)
6.  [Bottom line](#6-bottom-line)
7.  [Detailed plan for option B](#7-detailed-plan-for-option-b)
    1.  [Goals and non-goals](#71-goals-and-non-goals)
    2.  [Guiding principles](#72-guiding-principles)
    3.  [Package layout](#73-package-layout)
    4.  [Interface design](#74-interface-design)
    5.  [Target Python API (sketch)](#75-target-python-api-sketch)
    6.  [Work packages](#76-work-packages)
    7.  [Timeline](#77-timeline)
    8.  [Testing strategy](#78-testing-strategy)
    9.  [Parallelism and performance](#79-parallelism-and-performance)
    10. [Risks and mitigations](#710-risks-and-mitigations)
    11. [First 90 days](#711-first-90-days)
    12. [Decisions needed](#712-decisions-needed)

## 1. What the code is today

**Size and makeup** (Fortran line counts, not counting the bundled `finufft` and `h5fortran`)

| Component                    | Files | Lines     | What it does                                                                                                                                                                                                      |
|------------------------------|-------|-----------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `discus`                     | ~200  | **~156k** | Builds structures, simulates defects (MC/RMC/MMC, domains, stacking faults, molecules, surfaces, waves), and computes Fourier, powder, PDF, 3D-PDF and Patterson                                                  |
| `lib_f90`                    | ~150  | ~70k      | Shared base: the command interpreter, variables, macros, loops, expression evaluator, error lists, HDF5/NeXus I/O, FFT, random numbers, MPI/OpenMP glue, sockets, readline, and a Fortran↔Python bridge (`forpy`) |
| `kuplot`                     | 66    | ~38k      | Plotting (PGPLOT on X11, Windows or PNG) plus data loading, maths and least-squares fitting                                                                                                                       |
| `diffev`                     | 39    | ~10k      | Differential-evolution refinement, with MPI for parallel cost functions                                                                                                                                           |
| `refine`                     | 31    | ~8k       | Least-squares refinement that drives `discus` macros                                                                                                                                                              |
| `suite`                      | 16    | ~3k       | One umbrella program that joins all of the above and switches between their menus                                                                                                                                 |
| `mixscat`, `experi`, `nexus` | small | ~7k       | Extra tools                                                                                                                                                                                                       |

There is also C/C++: `finufft`, plus readline, sockets and glob helpers.

### Architecture

- **Everything runs through the command language.** Each program and each sub-menu has its own `*_loop` → `get_cmd` → `*_kdo`. In `discus`, 30 files contain such a menu loop, and the top-level `discus_kdo` dispatches about 70 commands. In `discus/prog` alone the code makes about 1,700 `str_comp` keyword comparisons, about 1,100 `get_params`/`ber_params` calls and about 3,400 `ier_num` error checks.
- **The command language is a whole small language.** It has variables (`i[]`, `r[]`, `res[]`), an expression evaluator, `do`/`if` blocks, macros with parameters, `branch` between programs, `system`, `socket` and `learn`.
- **All state is global.** Crystal, Fourier, powder, MC and similar data live in roughly 75 `*_mod` modules (for example `crystal_mod` holds `cr_pos`, `cr_iscat` and so on). Compute routines read their inputs from these modules and write results back into them, not through arguments.
- **Results come back as text and through `res_para`.** There are about 1,250 `WRITE(output_io…)` statements in `discus`, and numeric results go into `res[]` for macros to use.
- **The programs talk to each other through macros.** `refine` and `diffev` work by running user macros that call `discus` and `kuplot`, then reading the results back from variables or files. MPI spreads those macro runs across processes.
- **Graphics live in `kuplot`, which calls PGPLOT directly.** A no-graphics build already exists (`*_nographic.f90` and the `DIFFUSE_GRAPHICS=OFF` option).
- **Some useful pieces already exist:**
  - `python/suite_python.f90` (built with f2py) offers `initialize_suite`, `execute_macro`, `get_data` and `set_value`. Python can call it, but only by passing command strings.
  - `discus_py.f90`, `kuplot_py.f90` and `diffev_py.f90` expose `command()`, `get` and `send`.
  - `forpy` lets the Fortran side call Python.
  - There is a CUDA kernel for the structure factor, OpenMP in the Fourier, PDF and MMC code, and FINUFFT for non-uniform FFTs.
- **There are almost no tests.** The repo has only `sysmac/systest.mac` files, so there is no regression suite.

**What this means for a migration.** The compute kernels themselves are mostly decoupled. For example, `four_strucf` only uses `discus_config_mod`, `diffuse_mod` and `precision_mod`. The coupling sits in the menu layer around them: it parses commands, sets global module variables, calls the kernel, then prints the results or stores them in `res[]`. So the usual entry point is not "the algorithm" but "the global state that menu commands have built up". That global state is the main obstacle, more than the Fortran language itself.

## 2. Target architecture (shared by all options)

     Web UI / Jupyter / scripts / CLI
                │
       Python API layer   (Crystal, Structure, Fourier, Powder, PDF, MC, Fit objects)
                │                     │                    │
       I/O (CIF, HDF5/NeXus,     Plotting (matplotlib,   Refinement (scipy.optimize,
       ASE/pymatgen, xarray)     plotly, a web viewer)   lmfit, DE, dask/MPI)
                │
       Compute core: no global I/O and no parsing; arrays in, arrays out
       (the existing Fortran kernels behind a thin interface, or Python versions)

The rules this architecture follows:

- The core never prints, never parses text and never draws.
- The core reports errors through return codes or exceptions, not `ier_num` and a global error state.
- Plotting, file formats and refinement become ordinary Python libraries.
- The command language becomes optional: a compatibility layer at most, not the main interface.

## 3. High-level options

### A. Wrap the current interpreter (weeks)

Extend the existing `suite_python` and f2py bridge. Python sends command strings, and data comes back as NumPy arrays through `get_data`, `get_r` and `set_value`.

**Good:** almost no Fortran changes, works now, and is useful for Jupyter.

**Bad:** the command language remains the real API, global state stays, and you can't run two structures at once. Errors and output are still text. A web UI on top of this would be fragile.

**Best use:** a short-term stopgap, and a way to build a regression-test corpus.

### B. Fortran compute behind a clean Python API — recommended

1.  Find the computational entry points: structure factor and Fourier, powder, PDF, 3D-PDF, Patterson, MC/RMC/MMC moves and energies, symmetry and space-group expansion, stacking faults, domains, and so on.
2.  For each one, write a `bind(C)` wrapper that takes explicit arguments (positions, types, cell, q-grid, form factors) and returns arrays. At first, the wrapper can copy its arguments into the existing global modules and call the current routine. Later, refactor each routine so it takes arguments directly.
3.  Build it as a Python extension, either with f2py or with `bind(C)` + nanobind/pybind11/ctypes, packaged by scikit-build-core or meson-python. The second route handles derived types and allocatable arrays more cleanly and gives proper wheels.
4.  Write a Pythonic object model on top: `Crystal`, `Fourier(…)`, `Powder(…)` and similar, using NumPy and xarray.
5.  Replace `kuplot` plotting with matplotlib/plotly, `refine` with scipy/lmfit, and `diffev` with scipy's DE or a small Python DE that uses dask/mpi4py/Ray for parallel runs. Move file I/O to h5py, gemmi/ASE/pymatgen (CIF) and nexusformat.
6.  Optional: a thin "macro compatibility" translator for existing user scripts, or keep the old `discus_suite` binary as a legacy build alongside.

**Good:** keeps decades of validated physics and the existing speed (OpenMP, CUDA, FINUFFT). Each step can be tested against the old binary. It also leads naturally to a web or service backend.

**Bad:** removing the global state needs care. Several modules, especially `mmc` (~9k lines), `structur`, `chem` and `place_molecule`, are tightly tied to their menus. Building and packaging Fortran plus HDF5 and FFTW for wheels takes real work, though cibuildwheel can handle it.

### C. Full rewrite in Python (years)

Reimplement the physics in pure Python, with Numba, JAX or CuPy for the hot loops: structure factors, PDF histograms, MC.

**Good:** one language, simplest packaging, the easiest code for others to contribute to, and GPU support and automatic differentiation (JAX) for refinement.

**Bad:** about 160k lines of domain logic to re-derive and re-validate, with no test suite as a safety net. It would take years of effort and carries a high risk of subtle physics regressions.

**Best use:** only for selected parts, for example new features, or kernels that are short and easy to check, such as the Debye sum, PDF histogramming and the structure factor on a grid.

### D. Hybrid that grows into a rewrite (long term)

Start with B. Port each module to Python when it is cheaper to rewrite than to untangle: I/O, plotting, refinement, set-up and geometry tools, symmetry (using spglib/gemmi). Keep Fortran, C++ or CUDA for the heavy kernels (Fourier, MC energy loops, powder/PDF) for as long as they are faster.

This is option B followed by C, done as a long-term plan.

## 4. Suggested route (B, then D)

| Phase                             | Work                                                                                                                                                                                          | Outcome                                                  |
|-----------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------|
| **0. Safety net**                 | Use the existing `suite_python` (option A) to run a set of representative macros (Manual examples, user workflows) and save the outputs as golden HDF5/NumPy files. Add CI with pytest.       | Regression tests: the precondition for everything else   |
| **1. Packaging**                  | Add a `pyproject.toml` with scikit-build-core around the existing CMake. Build one importable `diffuse._core` with `DIFFUSE_GRAPHICS=OFF`, no X11 and no readline.                            | `pip install` works                                      |
| **2. Data model**                 | Python `Crystal`/`Structure` classes holding NumPy arrays, plus Fortran getters and setters for `crystal_mod` arrays (positions, types, cell, properties, molecules).                         | Structures move between Python and Fortran without files |
| **3. Kernel APIs**                | Fourier/structure factor, then powder (Debye and FFT), then PDF, then 3D-PDF/Patterson. Each gets a pure-argument entry point and is checked against the golden data.                         | The main calculations can be called from Python          |
| **4. Remove pieces from Fortran** | Replace plotting with matplotlib/plotly, refinement with lmfit/scipy plus Python DE, and file I/O with h5py/gemmi/ASE. Stop maintaining `kuplot`, `refine` and `diffev` as Fortran menus.     | Smaller Fortran code, Python-native workflows            |
| **5. Simulation engines**         | MC/MMC/RMC, domains, stacking, surfaces, molecules. These are the hardest because they hold the most state. Consider a callback or "move/energy plug-in" interface, or port parts to Numba.   | Defect modelling from Python                             |
| **6. Front ends**                 | Jupyter widgets first. Then a service layer (FastAPI with a job queue such as Celery, RQ or Dask, results in HDF5/Zarr) and a web UI (React with plotly or three.js for 3D reciprocal space). | Web front end                                            |
| **7. Legacy**                     | Freeze `discus_suite` as the classic interface, or offer a macro→Python translator for common commands.                                                                                       | Existing users keep working                              |

## 5. Main risks and decisions to make early

- **Global state vs. multiple instances.** A web service needs several independent sessions. With global Fortran modules, each session needs its own process, which is acceptable at first and naturally handled by workers. True re-entrancy requires refactoring into derived types.
- **Error handling.** `ier_num`/`ier_typ` and some `STOP` paths must become return codes turned into Python exceptions. A Fortran `STOP` would kill the Python process.
- **Output.** The many `WRITE(output_io…)` calls need a quiet or logging mode.
- **Precision and compiler matrix.** You currently support gfortran, ifort, PGI, OpenMP and optional CUDA. Wheels mean picking a supported subset.
- **User base.** Macros are how existing users work, so decide early whether you need backward compatibility (a translator, or the old binary frozen) or will make a clean break.
- **Licensing and dependencies.** PGPLOT disappears (a good thing). FFTW, HDF5 and FINUFFT stay as build dependencies.

## 6. Bottom line

The physics code is valuable and fairly modular, but it is wrapped in a command interpreter and global-state layer that makes up a large share of the code. The recommendation is **option B**, a Python API over the existing Fortran kernels, built in phases. Start by building the regression-test corpus and Python packaging, replace plotting, I/O and refinement with standard Python libraries, and port kernels to Python only where it pays off. That gives a scriptable Python package in months rather than years, and a natural backend for a web front end later.

## 7. Detailed plan for option B

Option B keeps the Fortran physics and replaces everything around it, one piece at a time. Each step must leave a working `discus_suite` and a passing test suite. The plan is written for a small team (about 1–2 developers plus a domain expert to review results). All durations are rough estimates.

### 7.1 Goals and non-goals

**Goals**

- `pip install diffuse` gives a Python package with no X11, PGPLOT or readline.
- Structures, Fourier, powder, PDF and defect simulations can be run from Python with NumPy arrays in and out. No command strings and no temporary files.
- Plotting, file I/O and refinement are handled by standard Python libraries.
- Errors become Python exceptions, and output goes to Python `logging`.
- The package works as the compute backend for a later web service.

**Non-goals (for option B)**

- Rewriting the numerical kernels in Python (that is option C/D, case by case).
- Making the Fortran core fully re-entrant in the first pass. Isolation comes from processes at first.
- Reproducing every command-language feature (loops, `branch`, `socket`, `learn`). Python replaces them.
- Dropping `discus_suite`. It stays as a legacy build until users have moved.

### 7.2 Guiding principles

- **Strangler pattern.** New Python entry points wrap existing routines. The command menus keep calling the same routines, so old and new paths share one implementation and both stay correct.
- **Golden data first.** No refactor lands without a regression test that compares against output from the current release.
- **Thin, flat C boundary.** The only thing that crosses from Fortran to Python is plain C types and contiguous arrays (`iso_c_binding`). No Fortran derived types, strings or allocatables cross it.
- **Python owns the data model.** The truth lives in Python objects. Fortran modules are a cache that gets filled before a computation runs. Later refactors can remove that cache one module at a time.
- **Keep the menu layer, stop growing it.** New features are added as Python API first, and as a menu command only if really needed.

### 7.3 Package layout

    DiffuseCode/
    ├── (existing Fortran trees: discus/ lib_f90/ kuplot/ diffev/ refine/ suite/ …)
    ├── capi/                        # NEW: Fortran bind(C) shims + C header
    │   ├── dc_capi_core.f90         #   init, reset, errors, log redirection
    │   ├── dc_capi_crystal.f90      #   get/set cell, atoms, scattering types, molecules
    │   ├── dc_capi_fourier.f90      #   single-crystal Fourier / structure factor
    │   ├── dc_capi_powder.f90
    │   ├── dc_capi_pdf.f90
    │   └── diffuse_capi.h
    ├── pyproject.toml               # NEW: scikit-build-core, drives the existing CMake
    └── python/diffuse/              # NEW: the Python package
        ├── _core.cpp                #   nanobind module over diffuse_capi.h
        ├── __init__.py
        ├── errors.py                #   DiffuseError hierarchy mapped from ier_num/ier_typ
        ├── crystal.py               #   Crystal / Structure / Molecule classes
        ├── io/                      #   CIF (gemmi), DISCUS .stru, HDF5/NeXus (h5py), xyz, …
        ├── scattering/              #   fourier.py, powder.py, pdf.py, patterson.py
        ├── defects/                 #   replace, domains, stacking, surface, waves, …
        ├── mc/                      #   MMC / RMC drivers
        ├── refine/                  #   lmfit + differential-evolution wrappers
        ├── plot/                    #   matplotlib / plotly helpers (replaces kuplot)
        └── legacy.py                #   run_macro() via suite_python, for migration

CMake already builds `lib_f90`, `discus_all` and the others as libraries, so the new extension can link them directly. One build step needs untangling: `discus` currently includes and links `kuplot`, so the dependency on `kuplot`'s plotting code has to be cut or replaced with stubs for the Python build. The existing `*_nographic.f90` files give a head start.

### 7.4 Interface design

#### The C-ABI shim

Each shim function takes plain arguments, fills the existing global modules, calls the existing routine and copies results back out. Example for the structure:

    module dc_capi_crystal
      use iso_c_binding
      use crystal_mod
      use discus_allocate_appl_mod
      implicit none
    contains
      integer(c_int) function dc_set_atoms(n, pos, iscat, prop, surf) &
           bind(C, name="dc_set_atoms")
        integer(c_int), value         :: n
        real(c_double), intent(in)    :: pos(3, n)     ! fractional coordinates
        integer(c_int), intent(in)    :: iscat(n), prop(n), surf(n)
        call dc_clear_error()
        if (n > NMAX) call alloc_crystal_nmax(n)       ! grow existing arrays
        cr_natoms        = n
        cr_pos(:, 1:n)   = pos
        cr_iscat(1:n,1)  = iscat
        cr_prop(1:n)     = prop
        ! … surface, molecule info, then update derived quantities
        dc_set_atoms = dc_status()                     ! 0 = ok, else ier_num
      end function
    end module

The same pattern applies to `dc_set_cell`, `dc_set_scat_types`, `dc_get_atoms` and others. The existing `discus_allocate_appl_mod` already contains the resizing routines (`alloc_crystal_nmax`, `alloc_crystal_scat`, `alloc_diffuse_four` and so on).

#### Calls into compute routines

- Each `dc_*_run` function sets the module variables that the menu would normally set (for example the corners, increments and radiation in `diffuse_mod`), then calls `four_run`, `powder_run`, `pdf_run` and so on.
- Results are copied into caller-supplied NumPy buffers (`csf`, `dsi`, powder and PDF arrays), or exposed through a getter that reports their shape.
- The menu `set` commands contain validation logic. Where it matters, move that logic into small shared routines so both the menu and the shim call them, instead of copying it.

#### Errors

- Every shim function starts by clearing `ier_num`/`ier_typ` and returns them at the end. The C++ layer turns non-zero values into a Python exception and takes the text from the existing error lists (`errlist`).
- About 30 `STOP` statements exist, for example in `four_strucf_mod.f90` and the `assert*` helpers in `utilities.f90`. Each one must become an error return. Otherwise a single bad input ends the Python process.
- Any path that would drop into an interactive prompt on error (as the current `discus_py` `command()` does) must be disabled when called from the API.

#### Output and logging

- Add a global "API mode" flag in `prompt_mod`. When it is set, `output_io` goes to a buffer or `/dev/null`, and selected messages are passed to a C callback that forwards them to Python's `logging`.
- This needs no changes to the ~1,250 `WRITE` statements; only the unit they write to changes.

#### State and sessions

- One Python process holds one Fortran state. `dc_reset()` calls the existing reset routines (`discus_reset_all`, `fourier_reset`, `powder_reset`, `pdf_reset`) so that each calculation starts clean.
- A Python `Session` object can run the core in a worker process (multiprocessing or loky). That gives isolation for notebooks running several models and for the web service, without making the Fortran re-entrant.
- Longer term, the existing `cl_cryst` class (`class_crystal.f90`) is a starting point for passing a crystal as a derived type instead of through `crystal_mod`.

### 7.5 Target Python API (sketch)

    import diffuse as dc

    cr = dc.io.read_cif("PbTiO3.cif")                 # or dc.io.read_stru("old.stru")
    sc = cr.expand(supercell=(20, 20, 20))            # space-group expansion + supercell

    rng = dc.Random(seed=42)
    sc = dc.defects.replace(sc, "Ti", "Zr", fraction=0.1, rng=rng)
    sc = dc.mc.mmc(sc, targets=[dc.mc.Correlation("Ti", "Zr", neighbour=1, value=-0.3)],
                   cycles=50, rng=rng)

    f = dc.scattering.single_crystal(sc, radiation="xray", wavelength=0.71,
                                     plane=dc.Plane(corner=(-5,-5,0), right=(5,-5,0),
                                                    top=(-5,5,0), npoints=(401,401)))
    f.intensity.plot()                                # xarray.DataArray with hkl coordinates

    p = dc.scattering.powder(sc, qmin=0.5, qmax=25, dq=0.01, method="debye")
    g = dc.scattering.pdf(sc, rmax=30, dr=0.01)

    fit = dc.refine.differential_evolution(model=my_model, data=obs,
                                           params={"frac": (0.0, 0.3)}, workers=8)

The results are `xarray` objects with coordinates and metadata, and they can be saved straight to HDF5/NeXus.

### 7.6 Work packages

### WP0 · Regression corpus

*Months 0–3 · prerequisite for everything else*

**Tasks:** Collect 30–60 representative macros (Manual examples, tutorials, real user workflows). Run them with the current release through `suite_python` and save structures, Fourier maps, powder patterns and PDFs as reference HDF5/NumPy files. Set fixed random seeds. Set up pytest and CI (GitHub Actions on Linux and macOS).

**Deliverable:** `tests/golden/` plus a CI job that runs it on every push.

**Done when:** All main command areas (read, symmetry, replace, domain, stack, fourier, powder, pdf, mmc, rmc) are covered by at least one test.

### WP1 · Build and packaging

*Months 1–4*

**Tasks:** Add `pyproject.toml` with scikit-build-core on top of the existing CMake. Add a `DIFFUSE_PYTHON_CORE` option that builds without X11, PGPLOT, readline or the `kuplot` plotting code. Build wheels with cibuildwheel for Linux (manylinux) and macOS, with FFTW and HDF5 bundled. Look at Windows support later.

**Deliverable:** `pip install .` and `import diffuse` both work, and CI builds wheels.

**Done when:** A clean virtual environment can install the wheel and run one smoke test.

### WP2 · Core infrastructure

*Months 2–6*

**Tasks:** Add the API-mode flag, output redirection and logging callback. Map errors to exceptions. Remove or convert every `STOP`. Add `dc_init` and `dc_reset`. Expose the random seed and state (`random_mod`) so runs are reproducible. Add a `legacy.run_macro()` helper for running old macros during the migration.

**Deliverable:** `diffuse._core` with init, reset, errors, logging and RNG control.

**Done when:** A deliberately bad input raises a Python exception. The process survives, prints nothing to the terminal, and the next call works.

### WP3 · Crystal data model and I/O

*Months 3–8*

**Tasks:** Python `Crystal` class (cell, space group, atoms, ADPs, occupancies, properties, molecules). Two-way transfer to and from `crystal_mod`. Readers and writers for CIF (gemmi), DISCUS `.stru`, xyz and HDF5. Space-group expansion through the existing `spcgr_*` routines, cross-checked with spglib.

**Deliverable:** `dc.Crystal`, `dc.io`.

**Done when:** Round-tripping all golden structures (Python → Fortran → Python) gives identical data, and the files written match the legacy `save` output.

### WP4 · Single-crystal Fourier (first vertical slice)

*Months 5–10*

**Tasks:** Shim over `four_run` and `four_strucf`: radiation, wavelength, corners and increments, form factors, ADPs, lots, Friedel symmetry and output arrays. Results returned as `xarray` with hkl coordinates. Optional builds with OpenMP, CUDA and FINUFFT.

**Deliverable:** `dc.scattering.single_crystal()`

**Done when:** Results match the golden data to a set tolerance (for example a relative difference of 1e-10 in double precision), and speed is within 5% of the binary.

### WP5 · Powder

*Months 8–13*

**Tasks:** Shims over `powder_run`, `powder_complete` and `powder_nufft` (Debye, complete, FFT). Profile functions, preferred orientation, instrumental settings. Output as Q, 2θ and d. Comparison with experimental data in Python.

**Deliverable:** `dc.scattering.powder()`

**Done when:** All powder golden tests pass.

### WP6 · PDF, 3D-PDF and Patterson

*Months 10–15*

**Tasks:** Shims over `pdf_setup`, `pdf_run` and `pdf_determine`, plus the 3D-PDF (`discus_3dpdf`), `exp2pdf` and Patterson routines. Keep the PDF-fitting code paths (`pdf_rmc_*`) for WP8.

**Deliverable:** `dc.scattering.pdf()`, `pdf3d()`, `patterson()`, `exp2pdf()`

**Done when:** All PDF golden tests pass.

### WP7 · Structure building and defects

*Months 9–18 · broad, so done in several releases*

**Tasks:** Python functions over the existing routines for symmetry operations, replace and remove, insert, domains (`domain`), stacking faults (`stack`), surfaces, waves, shear, molecules (molecularize, place molecules, decorate), and connectivity and neighbour analysis (`chem`, `conn`). Priority goes to the commands that users rely on most.

**Deliverable:** `dc.defects.*`, `dc.analysis.*`

**Done when:** Each wrapped command has a golden test. Stacking faults and domains are written up as notebook tutorials.

### WP8 · Monte Carlo (MMC/RMC)

*Months 14–22 · the most stateful part*

**Tasks:** Wrap `mmc` (~9k lines) and `rmc`. Express targets (correlations, distances, energies) as Python objects and send them to the existing set-up code. Report progress through a callback so the caller can follow and stop long runs. Later option: Python-defined energy terms through callbacks.

**Deliverable:** `dc.mc.mmc()`, `dc.mc.rmc()`

**Done when:** With fixed seeds, results are statistically the same as the golden runs, and runs can be stopped and resumed.

### WP9 · Refinement and plotting replaced in Python

*Months 6–14 · can run alongside WP4–WP7*

**Tasks:** Least-squares refinement with lmfit (replaces `refine`). Differential evolution (replaces `diffev`), either scipy's or a ported DE that keeps the current DIFFEV options, with parallel workers through concurrent.futures, mpi4py or Dask. Plot helpers in matplotlib and plotly for 1D, 2D maps and 3D slices (replaces `kuplot` for typical uses). Readers for the formats `kuplot` loads.

**Deliverable:** `dc.refine`, `dc.plot`

**Done when:** Two or three published refinement examples are reproduced in notebooks.

### WP10 · Service layer and web front end (first version)

*Months 16–24*

**Tasks:** A FastAPI service with a job queue (Celery, RQ or Dask) and one worker process per job. Pydantic schemas for model and calculation requests. Results stored as HDF5/Zarr. A first web UI: load a structure, set up a single-crystal or powder calculation, view the results interactively (plotly; three.js for 3D).

**Deliverable:** A deployable container image and a demo UI.

**Done when:** A user can go from a CIF file to a diffuse-scattering map in a browser.

### 7.7 Timeline

Months from the start, assuming about 1.5 full-time developers. Some packages overlap because they can be worked on in parallel.

| Work package                       | Months 0 → 24 (one block = one month)                                                         |
|------------------------------------|-----------------------------------------------------------------------|
| WP0 Regression corpus              | `███·····················` 0–3      |
| WP1 Build and packaging            | `·███····················` 1–4    |
| WP2 Core infrastructure            | `··████··················` 2–6    |
| WP3 Crystal model and I/O          | `···█████················` 3–8   |
| WP4 Single-crystal Fourier         | `·····█████··············` 5–10  |
| WP9 Refinement and plotting        | `······████████··········` 6–14    |
| WP5 Powder                         | `········█████···········` 8–13  |
| WP7 Structure building and defects | `·········█████████······` 9–18  |
| WP6 PDF, 3D-PDF, Patterson         | `··········█████·········` 10–15 |
| WP8 Monte Carlo                    | `··············████████··` 14–22 |
| WP10 Service and web UI            | `················████████` 16–24 |

**Milestones:**

- **M1 (month 6):** installable package, error handling and logging, crystal round trip. Alpha release for internal users.
- **M2 (month 12):** single-crystal Fourier, powder, Python refinement and plotting. First public beta, usable for most diffuse-scattering workflows from notebooks.
- **M3 (month 18):** PDF and the main defect tools. Most of what users do no longer needs macros.
- **M4 (month 24):** Monte Carlo and the first web front end. Python becomes the primary interface; `discus_suite` moves to maintenance only.

### 7.8 Testing strategy

- **Three levels.** (1) Golden tests against the legacy binary. (2) Unit tests for the shims (shapes, errors, reset). (3) Physics checks that don't depend on legacy output: Bragg positions, Friedel symmetry, the Debye result vs. the Fourier result for the same structure, PDF normalisation, sum rules.
- **Tolerances.** Use tight tolerances for deterministic calculations. Use statistical comparisons (means and spreads over seeds) for MC/RMC.
- **CI matrix.** gfortran on Linux and macOS for every push; ifort/ifx and the CUDA build nightly or on demand.
- **Performance tracking.** Benchmark key kernels (asv or pytest-benchmark) so the wrapper layer never adds noticeable overhead.

### 7.9 Parallelism and performance

- **OpenMP** continues to run inside the kernels. The nanobind layer releases the Python GIL during long calls, so notebooks stay responsive and threads can run other work.
- **MPI**, used today by `diffev` and `refine` to run macros in parallel, is replaced by Python process pools (concurrent.futures, Dask, Ray or mpi4py on clusters). Each worker holds its own Fortran state, which fits the global-module design.
- **CUDA and FINUFFT** become optional build variants, for example `diffuse[cuda]`. Python picks the backend at run time and falls back to CPU.
- **Data copies.** Arrays cross the boundary as contiguous NumPy arrays (Fortran order where needed). Large outputs such as 3D Fourier volumes are written straight into buffers that Python allocates, so there is no extra copy.

### 7.10 Risks and mitigations

| Risk                                                                      | Impact                                                 | Mitigation                                                                                                                                                                   |
|---------------------------------------------------------------------------|--------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Hidden coupling: kernels depend on module state set by unrelated commands | Wrong results that look plausible                      | Golden tests per command. Start each shim call from `dc_reset()` and set every relevant variable explicitly. Write down which module variables each kernel reads.            |
| Validation logic lives only in menu `set` handlers                        | The API accepts invalid input                          | Move validation into shared routines that both the menu and the shim call. Add pydantic validation in Python.                                                                |
| Packaging Fortran with HDF5 and FFTW as wheels                            | Installation problems, especially on macOS and Windows | cibuildwheel with pinned dependencies. Offer a conda-forge package as well, which is often easier for scientific users.                                                      |
| Development continues in the Fortran code during the migration            | Merge conflicts, a moving target                       | Keep shims in a separate `capi/` tree. Golden tests catch behaviour changes. Agree that new features are added through the API first.                                        |
| Users depend on macros                                                    | Slow adoption                                          | Keep `discus_suite`. Provide `legacy.run_macro()`. Write side-by-side "macro → Python" guides, and possibly a translator for common commands.                                |
| Size of MMC/RMC                                                           | WP8 overruns                                           | Schedule it last, wrap rather than refactor, and limit the first release to the most-used energy terms.                                                                      |
| Limited developer time and few people with Fortran experience             | The plan slips                                         | Release in milestones so value arrives early. Make the Python layer welcoming to new contributors. The Fortran work is concentrated in the shim and infrastructure packages. |

### 7.11 First 90 days

- [ ] Agree on package name, licence, supported platforms and minimum Python version.
- [ ] Choose and document 30+ golden macros. Run them with the current release and store the outputs.
- [ ] Set up pytest and GitHub Actions, run the golden tests through `suite_python`.
- [ ] Add the `DIFFUSE_PYTHON_CORE` CMake option and cut the `discus` → `kuplot` plotting dependency for that build.
- [ ] Add `pyproject.toml` (scikit-build-core) and a minimal nanobind module that only calls `dc_init` and returns the version.
- [ ] Prototype the API-mode output redirect and error-to-exception mapping.
- [ ] Make an inventory of all `STOP` statements and interactive fallback paths, and a plan to remove each one.
- [ ] Prototype `dc_set_cell`, `dc_set_atoms` and `dc_get_atoms`, plus a Python `Crystal` that round-trips a golden structure.
- [ ] Write down which module variables `four_run` reads, as preparation for WP4.

### 7.12 Decisions needed

- **Binding technology:** `bind(C)` + nanobind (recommended: cleaner arrays, exceptions and GIL control) or f2py (already used by `suite_python`, less code at first, weaker for complex interfaces).
- **Distribution:** PyPI wheels, conda-forge, or both.
- **Supported platforms and compilers** for official builds. Is Windows native, or WSL only?
- **Backward-compatibility promise** for macros and file formats, and how long `discus_suite` is maintained.
- **Data formats:** NeXus/HDF5 as the standard result format, and which version of the DISCUS structure format to keep.
- **Web goals:** local single-user app, shared group server or public service. This decides how much work goes into authentication, quotas and job management in WP10.

Based on a read-only review of the DiffuseCode repository (version 6.23.01). Line counts are approximate.

