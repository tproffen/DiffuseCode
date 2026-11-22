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

### Keyboard input: 18 reads from standard input (22 on a closer look, see step 2)

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

1. **Output.** *Done (Nov 2026), differently from first planned.* A log
   file opened by the Fortran code would be buffered by gfortran (text
   arrives in blocks of several kB, so no live progress), and the only
   switch against that (`GFORTRAN_UNBUFFERED_ALL`) would slow down every
   large file DISCUS writes. Standard output, by contrast, is unbuffered
   (`GFORTRAN_UNBUFFERED_PRECONNECTED`, set by `_session.py`). So the
   file-descriptor capture stays as the output channel, and:
   * the banner, the start-up text and the `cd` echo write to `output_io`
     (standard output in `discus_suite`, so unchanged there); in library
     mode an error writes all its lines (also DISCUS's extra lines) to
     `error_io`, without the copy to the screen;
   * in library mode `set prompt` cannot move the output to a log file or
     switch it off;
   * `execute_macro` releases the Python GIL (`!f2py threadsafe`), so a
     thread can follow the output: with `echo=True` lines are shown as
     DISCUS writes them, in a terminal and in a Jupyter cell (tested with
     a kernel: lines 1 s apart arrive 1 s apart).
   The remaining `write(*` (debug dumps, rare internal messages) are still
   caught by the capture.
2. **No keyboard reads in library mode.** *Done (Nov 2026).* A closer
   look found 22 reads from standard input outside comments: 3 in the
   separate socket client `remote.f90` and 4 in main programs (not part
   of the library); the command reader (`get_cmd`, 2) is already safe
   (from Python, `linteractive` is false and commands come from the
   macro).
   In library mode:
   * `wait` and `wait return` continue at once; `wait input` is a new
     error, -21 COMM "Command needs a terminal (library mode)";
   * the pager pauses (`help`, `show`) do not wait (`help` and `manual`
     are switched off without a terminal anyway); the "remainder of macro
     might be faulty" prompt answers "cancel";
   * `exit` at the top level ends the macro only: no program exit, so no
     clean-up and no "Update Ubuntu now?" question, which treated end of
     input as yes and would have run the update;
   * as a safety net `_session.py` makes standard input, during each call,
     a file of empty lines: any other read gets an empty line, as if
     Return was pressed (the `exp2pdf` questions keep their defaults).
     Not `/dev/null`: after an end of file gfortran stops the whole
     process at the next read of that unit (tested).
3. **Commands without a macro file.** *Done (Nov 2026).* New library
   routine `macro_store` (`lib_f90/macro_internal.f90`) puts lines into
   the macro storage under a name; `dc_run_commands(text)` stores the
   commands as `/diffuse/python/commands.mac` and runs them through the
   usual macro machinery (`file_kdo` looks into the storage before the
   disk), so loops, blocks, nested `@file` macros, error lines and
   `macro_close` work as before. `run()` writes no file; older module
   builds still use `_driver.mac`. Lines longer than 1024 characters
   (`PREC_STRING`) are refused in Python.
4. **Reset completeness.** *Done (Nov 2026).* The `04.sro3` case was not
   left-over state but a bug in the bundled finufft (a chunk of the kernel
   coefficients stayed uninitialised when OpenMP ran fewer threads than
   planned, e.g. `OMP_THREAD_LIMIT=1`); fixed, see the WP0 log. Test:
   after an MMC run and `reset()`, a Fourier map and a powder pattern are
   identical to a fresh process (serial; with several threads the powder
   sums differ in the last digits from run to run anyway). The `diffuse`
   functions set every option they use, so they do not depend on `reset`.
   `reset` keeps variables, the random generator and `set` options such as
   `set parallel` (a process-wide `OMP_SET_NUM_THREADS`) and `set error`.
   Also found: `random_state()` missed the Gaussian number `gasdev` keeps
   waiting (`gset`, a SAVEd local); it is now in `random_mod`, and the
   state restores exactly also between the two numbers of a pair.
5. **Decide the binding** (plan: f2py now, `bind(C)` + nanobind or ctypes
   later). The `dc_*` functions are written so that this is a mechanical
   change; f2py stays until WP3/WP4 need arrays of structures.
