# WP2 – Core infrastructure: survey and first step

Survey of the Fortran code for WP2 (migration plan, section 7), done in
Nov 2026 on `claude-python-app`, and what the first step changed. Counts
cover `discus`, `kuplot`, `suite`, `diffev`, `lib_f90`, `refine` and
`experi`, without comment lines and without the bundled `finufft`,
`h5fortran` and `forpy_mod.F90`.

WP2 goal (plan): *a deliberately bad input raises a Python exception; the
process survives, prints nothing to the terminal, and the next call works.*
With the first step below this holds for command errors (test:
`tests/python/test_wp2_core.py`). Two gaps remain: "prints nothing" still
relies on `_session.py` capturing the output file descriptors, and the few
commands that read from the keyboard would wait for it (steps 1 and 2 in
section 3).

## 1. Findings

### STOP statements: 31, few reachable from Python

| Where | Count | Reachable from Python | Status |
|---|---|---|---|
| `errlist` when the error status is `exit` (`set error,exit`) | 1 | yes: one command ended the Python process | **fixed**: no STOP in library mode |
| menu loops (`suite`, `discus`, `kuplot`, `diffev`, `refine`, `experi`): read error on standard input (error -9, I/O) | 6 | yes, when stdin closes or fails | **fixed**: an error like any other in library mode |
| `four_strucf_version_15`: allocation failure | 1 | yes, for very large Fourier grids | **fixed**: error -3 (COMM), "Could not allocate arrays" |
| main programs (`discus.f90`, `kuplot_exit.f90`, `diffev.f90`, `suite.f90` 4x) | 7 | no (not part of the library) | leave |
| `utilities.f90` assert/nrerror (Numerical Recipes helpers) | 10 | only on a programming error (inconsistent array sizes), not on bad input | leave for now; convert when those routines are touched |
| `kuplot/gsas.f90` | 4 | GSAS file reading | convert with WP9 (KUPLOT) |
| `diffev/with_mpi.f90` | 1 | MPI only | leave |
| `lib_f90/fftpack6.0.f90` | 1 | library internal error | leave |

### Errors: one good path already

* 1054 places raise an error by setting `ier_num` (negative) and `ier_typ`
  and returning; the menu loops then call `errlist`, which prints the
  message through `disp_error`. That single exit point is where the first
  step records the error for Python.
* 88 places call `errlist` themselves, in the middle of a calculation
  (mostly warnings). These print but do not stop the command; Python sees
  them only in the output text, unless they are errors (then they are
  recorded too, but a command that carries on does not raise).
* Error classes (`errlist_mod`): COMM, FORT (expressions), I/O, MAC, MATH,
  APPL and per-program ones (CHEM, FOUR, MMC, PDF, RMC, HDF).

### Output: almost all through one unit

* 1668 writes go to `output_io`, 12 to `error_io`: one variable each, so a
  later step can send all of it to a file or a Python callback instead of
  capturing file descriptors.
* 357 writes go to `*` (standard output directly): 133 in
  `unified_dump.f90` (debug dumps), the rest in the start-up banner
  (`suite.f90`), `appl_unix.f90`, `exit_mod.f90`, KUPLOT's NeXus reader,
  and `disp_error` itself, which writes every error twice (to `error_io`
  and to `*`). To be moved to `output_io` before the descriptor capture in
  `_session.py` can go.

### Keyboard input: 18 reads from standard input

`wait`, `help`, parts of `show`, `remote`, `exit_mod`, `exp2pdf` and the
command reader `sup.f90` (when no macro is active). In library mode these
must not wait for the keyboard: to be turned into errors (or no-ops) in a
later step. Through the Python API they are only reached by `run_macro`
with such commands.

### Global state

About 460 modules hold state at module level (discus 192, lib_f90 150,
kuplot 61, diffev 32, refine 27). The suite command `reset` calls the
existing reset routines of DISCUS (`discus_reset_all`, 26 menus), KUPLOT,
DIFFEV and REFINE; it keeps variables (`r[]`, `i[]`, user variables) and
the random generator. The `04.sro3` failure (WP0 log: broken Fourier after
an MMC run, from state left in memory) shows that some state survives
where it should not; finding it is part of WP2 / the DISCUS bug batch.

### Random numbers

One generator, the compiler's `RANDOM_NUMBER`, with the state in
`random_state_mod` (`seed_vals`) plus the cached second Gaussian in
`random_mod` (`iset`). `ini_ran_ix` seeds it (the `seed` command; since
WP0 every word of the state is set from the user's seed). With several
OpenMP threads the order of draws can differ between runs, so exact
repeats need one thread (the tests use `OMP_THREAD_LIMIT=1`).

### How Python drives DISCUS today

`_session.py` writes the commands to a driver macro file and runs it with
`execute_macro`, which reads it through the macro machinery (macros are
kept in memory by name and emptied after each top-level run). A completion
flag (`r[499]`) on the last line tells whether a command failed. This
stays the transport for now: menus (read, fourier, output ...) read their
commands line by line from the macro, so a direct "run this one command"
call needs a line source for the sub-menus first.

## 2. First step (this branch)

Fortran (`lib_f90/errlist_mod.f90`, `errlist.f90`, the menu loops,
`four_strucf_mod.f90`):

* `ier_api` (library mode): no STOP on errors (cases above marked fixed).
  Off by default, so `discus_suite` behaves as before.
* `ier_last_*`: `disp_error` keeps the first error since it was cleared:
  number, type, type name, standard message and DISCUS's extra lines.

Python interface (`python/suite_python.f90`), functions with plain
integer, string and integer-array arguments only, so they can move to a
`bind(C)` layer later without changing the Python side:

| Function | Does |
|---|---|
| `dc_init` | `initialize_suite` in library mode |
| `dc_reset` | the suite command `reset` |
| `dc_clear_error`, `dc_last_error` | the recorded error |
| `dc_seed` | the command `seed n` |
| `dc_random_size`, `dc_get_random_state`, `dc_set_random_state` | save and restore the generator state |

Python (`diffuse`): `MacroError` has `code`, `category`, `text` and
`details`, taken from DISCUS's record instead of parsing the output;
`Session.reset`, `seed`, `random_state`, `set_random_state`,
`last_error`; `diffuse.seed()` and `diffuse.reset()`. Builds of
`suite_python` without the `dc_*` functions still work (old behaviour).

## 3. Next steps for WP2, in order

1. **Output to a callback.** Move the `write(*` output of the banner,
   `disp_error` (the duplicate) and `exit_mod` to `output_io`; in library
   mode open `output_io` on a file that `_session.py` reads, or call a
   Python logging function. Then drop the file-descriptor capture.
2. **No keyboard reads in library mode** (the 18 above): an error saying
   the command needs a terminal.
3. **Commands without a macro file.** Feed the command lines to the macro
   machinery from memory (it already stores macros in memory), so
   `run()` no longer writes `_driver.mac`.
4. **Reset completeness.** Find what `reset` leaves behind (start with the
   `04.sro3` MMC case, debug build with `-fcheck=all`).
5. **Decide the binding** (plan: f2py now, `bind(C)` + nanobind or ctypes
   later). The `dc_*` functions are written so that this is a mechanical
   change; f2py stays until WP3/WP4 need arrays of structures.
