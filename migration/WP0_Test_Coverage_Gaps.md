# Test coverage gaps – what to add and harden

Part of WP0 ([`WP0_Test_Baseline.md`](WP0_Test_Baseline.md)) of the [migration plan](DiffuseCode_Python_Migration_Plan.md). This compares the commands
documented in the help files with what the current macro tests actually run
and check. It then lists the tests to add or strengthen, ordered by when the
migration work packages (WP3–WP9) need them. `mixscat` is deprecated and not
included.

Generated from the run of 2026-11-13 (`claude-python-app` after PR #33; 31
tests). Re-run at any time:

```bash
python3 tools/run_legacy_tests.py --exe /path/to/discus_suite -j 4 --workdir /tmp/run
```

```bash
python3 tools/command_coverage.py --logs /tmp/run --markdown coverage.md
```

---

## 1. Summary

| Program | Commands reached | Sub-commands / options reached | Comment |
|---|---|---|---|
| DISCUS | 28 / 61 | 136 / 402 | Core scattering paths covered; I/O, symmetry, many structure-editing commands, RMC and Patterson not reached |
| KUPLOT | 34 / 82 | 9 / 69 | Used mainly as a test harness (`load`, `rvalue`) and for plot styling; `fit` and data processing not reached |
| DIFFEV | 0 / 28 | – | **Not tested at all** |
| REFINE | 6 / 9 | 3 / 4 | One test |
| SUITE | 3 / 7 | – | `diffev` branch, `parallel` not reached |
| Command language (`lib_f90`) | 19 / 31 | – | Covered as a side effect; `seed` never used |
| EXPERI | 0 / 1 | 0 / 10 | `exp2pdf` not reached (still in use?) |

**Reached is not the same as checked.** Many tests run a command but only
check something simple afterwards. Section 3 rates each test's check.

How the numbers are made, and their limits:

* Commands come from the help files: level 2 is a command, level 3 is a
  sub-command or option.
* What ran comes from the command echo in the test logs, matched by unique
  prefix, since the command language accepts abbreviations.
* Macro echo also prints the lines of `if` branches that are not taken.
  Commands inside `do` loops in a sub-menu are echoed with the outer prompt.
  A few counts can therefore be off; for example, `transform` counts as
  reached because a `trans …` line inside a menu loop matched it.
* Options passed as parameters (for example `set technique:nufft`) are only
  counted for `set`. Everything else is counted at command level.

---

## 2. Gaps by migration work package

### WP3 – Crystal model and I/O (needed first)

| Gap | Today | Why it matters |
|---|---|---|
| Structure export/import: `export` (SHELX, CIF, RMCProfile, VASP …), `import` | Not reached. `01.Simulation/export.mac` exists with 9 checks but is **not** in `test_all.mac`, and it fails 6 of 9 (see 4.1) | The Python `Crystal` class must read and write the same files |
| `save` options (`format`, `incl`, `omit`, `prop`, `asym` …) | 4 / 13 options | Round trip Python → Fortran → file |
| `spacegroup`, `wyckoff`, `asym`, `symmetry` details | Not reached / 14 of 24 sub-commands | Space-group expansion is wrapped in WP3 |
| `transform` (new cell, origin shift) | No sub-commands reached | Cell transformations are part of the data model |
| `d2r`, `r2d`, `index2cell`, `vprod` | Mostly not reached | Small metric helpers; cheap exact tests |
| KUPLOT `load` formats (`csv`, `h5`, `nexus`, `mrc`, …), `ksave`, `dsave` | 3 / 16 formats | Data-file parity when I/O moves to h5py |

### WP4 – Single-crystal Fourier

| Gap | Today | Why it matters |
|---|---|---|
| Radiation and energy: `electron`, `energy` (instead of `wvle`) | Not reached | Every radiation type needs a golden result |
| Anomalous dispersion `delf`, single `hkl`, `layer` | Not reached | Options the API must expose |
| `set technique` (standard vs. NUFFT) | Used once | The two methods should agree: a cheap consistency test |
| Only `set aver` and `set technique` used | – | Check the other `set` options (filters, lots variants, ADPs) |
| `diff-four` | Not reached | |
| Strong checks exist (R < 0.001 against `RESULTS/`) | 02.1d, 02.1d_periodic, 02.2d | Good model for the rest |

### WP5 – Powder

| Gap | Today | Why it matters |
|---|---|---|
| Neutron and electron powder (`neut`, `electron`) | Not reached | |
| One parameter set only | All 3 powder runs use the same `set` values (axis, calc, profile, lpcor, …) | Need variations: Debye vs. complete vs. FFT method, profile shapes, preferred orientation, 2θ vs. Q output |

### WP6 – PDF, 3D-PDF, Patterson

| Gap | Today | Why it matters |
|---|---|---|
| **`03.pdf` has no check** | Runs 9 PDF calculations; nothing is compared | Highest priority: add a G(r) comparison against a reference |
| PDF fitting mode (`pdf run`, `data`, `i/jsele`) | Not reached | |
| `patterson` | Not reached | Wrapped in WP6 |
| 3D-PDF | `RESULTS/*.3DPDF` exist in SRO_2, SRO_3, SRO_4, SRO_MOLECULES but are not compared | Easy win: compare them in the SRO tests |
| `exp2pdf` (DISCUS) and EXPERI `exp2pdf` | Not reached | Clarify whether EXPERI is still used |

### WP7 – Structure building and defects

| Gap | Today | Why it matters |
|---|---|---|
| Atom editing: `remove`, `kick`, `displacement`, `change`, `copy`, `define`, `append`, `anis` | Not reached | Small deterministic tests with exact atom-position checks |
| `insert` options (`fuzzy`, `density`, `origin`, axes …) | Only the plain form | |
| Molecules: `molecularize`, `demolecularize`, molecule plotting | Not reached; `04.SRO/SRO_MOLECULES` exists but is not in the suite and aborts (see 4.1) | |
| `perioditize`, `superspace`, `storage`, `switch`, `project`, `inverse`, `bragg` | Not reached | |
| `stack` modes `random`, `four`, `rotate`, `sigma` | Not reached | |
| `waves` (`osci`, `rot`, `shift`, `dens` …) | 15 / 26 | |
| `chemistry` analysis (`neig`, `env`, `homo`, `angle`, …) | 7 / 18 | Analysis output is exact and easy to compare |
| `connectivity` | Reached, no sub-commands | Used by decoration and molecules |

### WP8 – Monte Carlo

| Gap | Today | Why it matters |
|---|---|---|
| `rmc`, `mc` (old MC) | Not reached | Wrapped in WP8 |
| `mmc` | 5 tests, all with statistical tolerances (R < 0.25–0.6, correlation ±0.04–0.09) | Needs `seed` fixed (WP0 task 0.5a); then exact golden results |
| `seed` | Never used in a test | Fix first (0.5a), then use it everywhere |

### WP9 – Refinement and plotting

| Gap | Today | Why it matters |
|---|---|---|
| **DIFFEV: 0 / 28 commands** | Not tested. Candidates exist in `07.Refinement/EXAMPLE1`, `EXAMPLE_COMP`, `arctan` | The Python DE replacement needs reference behaviour |
| REFINE `fix`, `sigma`, `finished` | Not reached | |
| KUPLOT `fit` (least squares) | Not reached | Compare against lmfit/scipy in WP9 |
| KUPLOT data processing: `smooth`, `rebin`, `derivative`, `integral`, `fft`, `merge`, `convolute`, `mean`, `spline`, `sort` | Mostly not reached | **Decide:** if these move to NumPy, test only what users rely on |
| KUPLOT plot styling (colours, fonts, frames, titles) | Reached but not checked | **Do not test.** Replaced by matplotlib/plotly |
| MPI runs (`09.Nanoparticle_refinement`) | Disabled | Nightly or manual job |

### Command language and suite

The command language is being retired, so test only what macros still depend
on. `if`, `do`, `variable`, `fopen`, `fput` and similar are already covered as
a side effect. Exceptions:

* `seed`: needed for reproducible tests (0.5a).
* Error handling: no test checks that a bad command gives the expected
  error. WP2 maps errors to Python exceptions, so a small set of
  "this input must fail with error N" tests is worth having.

---

## 3. Strength of the existing checks

| Check | Tests | Assessment |
|---|---|---|
| **Exact reference data**, tight tolerance (R < 0.001–0.01 against `RESULTS/`) | 02.1d, 02.1d_periodic, 02.2d, 02.powder, 08.nano_simple, 08.nano_ellipsoid, 12.waves_density, 13.analyze | Strong. Model for new tests |
| **Statistical / loose** | 02.2d_lots (R < 0.13–0.15), 04.sro0/1/3/4 (R < 0.25–0.6), 05.stack_growth, 06.perovskite, 07.refine, 06.triangle (intensity ratios) | Catches large breaks only. Tighten after `seed` works |
| **Structural only** (atom counts, a few positions, cell size) | 01.circle, 01.checkerboard, 01.symmetry, 01.shear, 05.stack, 05.stack_file, 05.stack_types, 06.cube, 06.tilt, 11.properties, 12.waves (internal consistency), 14.deco01_all, 14.deco02 | Checks the structure is plausible, not that it is right. Add golden structure files (WP0 task 0.6) |
| **None** | 03.pdf | Add a check (task 0.4) |

Plot macros now write CIF files (`program cif`), but nothing compares them.
`08.Nanoparticle_example/SIMPLE/PLOT/block.cif` came out byte-identical to the
committed copy, so these could become cheap golden checks.

---

## 4. Work list

### 4.1 Repair existing tests that are not in the suite

These macros check themselves, but `test_all.mac` does not run them. Run
through the new runner, none of them passes today:

| Macro | Result | Likely cause |
|---|---|---|
| `01.Simulation/export.mac` (9 checks) | 6 FAIL, 3 PASS | SHELX and CIF exports differ from `Original/` in 2 atoms. Either the export changed (**possible regression, needs a human look**) or the reference files are outdated. The macro also writes `cmm2_export.rmc6f` but compares `cmm2_export.rmcf6` (typo) |
| `04.SRO/SRO_2/sro.2.mac 0.7,-0.7,1.5` | FAIL | Commented out in `test_all.mac`; the result is outside its tolerance |
| `04.SRO/SRO_MOLECULES/sro.2.mac 0.7,-0.7,1.5` | ABORTED: error opening a file | Input file missing or renamed |
| `08.Nanoparticle_example/Q/nano.mac 30, zno` | ABORTED: KUPLOT "Incompatible data sets for KCAL" | Output format changed since the reference was made |

### 4.2 Harden existing tests

| # | Task | WP0 task |
|---|---|---|
| H1 | `03.pdf`: compare G(r) against a stored reference | 0.4 |
| H2 | SRO tests: also compare the existing `RESULTS/*.3DPDF` files | 0.4 / 0.6 |
| H3 | Store golden outputs for the "structural only" tests (final structure, CIF) and compare atom positions | 0.6 |
| H4 | After `seed` is fixed, run MC tests with a fixed seed and tighten tolerances (SRO, stack growth, domains, 2d_lots, refine) | 0.5a / 0.5 |
| H5 | Compare the CIF files written by the plot macros | 0.6 |
| H6 | Run deterministic tests with 1 and N OpenMP threads and compare | 0.5 |

### 4.3 New tests, by priority

Ordered by when the migration needs them. Each test should be small (under
10 s), deterministic (fixed seed), compare one feature against a golden result,
and need no interactive viewer.

| ID | Area | Test | Check | Needed by |
|---|---|---|---|---|
| A1 | I/O | `read` → `save` → `read` round trip for `.stru`/`.cell`/free cell; all `save` options | exact atoms | WP3 |
| A2 | I/O | `export`/`import` CIF, SHELX, RMCProfile, VASP (repair `export.mac` first) | file + atoms | WP3 |
| A3 | Symmetry | `spacegroup` expansion for representative groups (triclinic … cubic, non-standard settings), `wyckoff`, `asym` | atom count + positions | WP3 |
| A4 | Metric | `transform` to a new cell and back; `d2r`, `r2d`, `vprod`, `index2cell` | exact numbers | WP3 |
| A5 | Fourier | X-ray / neutron / electron, `energy` vs `wvle`, `delf`, single `hkl` | golden intensity | WP4 |
| A6 | Fourier | Standard vs NUFFT technique on the same structure | agree within tolerance | WP4 |
| A7 | Fourier | 3D volume and lots variants | golden | WP4 |
| A8 | Powder | Neutron and electron; Debye vs complete vs FFT; profiles; preferred orientation; 2θ vs Q | golden pattern | WP5 |
| A9 | PDF | X-ray and neutron G(r) with damping/broadening options | golden G(r) | WP6 |
| A10 | Patterson / 3D-PDF | `patterson`; 3D-PDF of a small disordered crystal | golden map | WP6 |
| A11 | exp2pdf | DISCUS `exp2pdf` on a stored powder pattern | golden G(r) | WP6 |
| B1 | Editing | `remove`, `kick`, `displacement`, `change`, `copy`, `define`, `append`, `anis` | exact atoms | WP7 |
| B2 | Insert | `insert` options (origin, axes, density, fuzzy) | exact atoms | WP7 |
| B3 | Molecules | `molecularize`/`demolecularize`; repair SRO_MOLECULES | molecule table | WP7 |
| B4 | Defects | `perioditize`, `superspace`, `storage`, `switch`, `project`, `inverse`, `bragg` | exact / golden | WP7 |
| B5 | Stack | `random`, `four`, `rotate`, `sigma` modes | golden | WP7 |
| B6 | Analysis | `chemistry` `neig`, `env`, `homo`, `angle`, `bval` on a fixed structure | exact numbers | WP7 |
| B7 | Waves | remaining wave types (`osci`, `rot`, `shift`, `dens`) | golden | WP7 |
| C1 | MC | `mmc` with fixed seed: exact correlations (after 0.5a) | exact | WP8 |
| C2 | MC | `rmc` small fit to a stored pattern | R-value + seed | WP8 |
| C3 | MC | old `mc` (only if still in use) | | WP8 |
| D1 | DIFFEV | small DE fit (polynomial or `arctan` example) with seed, serial | final parameters | WP9 |
| D2 | REFINE | `fix`, `sigma`, `finished` | parameters + errors | WP9 |
| D3 | KUPLOT | `fit` of a known function | parameters | WP9 |
| D4 | KUPLOT | data processing actually used by users (after the keep/drop decision) | golden arrays | WP9 |
| E1 | Errors | Bad commands, bad parameters and missing files give the expected error number; the run continues | error code | WP2 |
| E2 | Command language | `seed` reproducibility: same seed → same random numbers | exact | 0.5a |

### 4.4 Out of scope (decide once, then document)

* **KUPLOT plot styling** (colours, fonts, frames, titles, marks): replaced by
  matplotlib/plotly, so not worth testing.
* **Command-language features** not used by any macro (`learn`, `mount`,
  `fsub`, `timer`, sockets): retire rather than test.
* **mixscat**: deprecated.
* **EXPERI**: test only if it is still used.

---

## 5. Conventions for new tests

* One feature per macro, small crystal, run time under 10 s.
* `seed` at the start of every stochastic test (once 0.5a is fixed).
* Compare against stored reference data (`RESULTS/` now, `golden/` from task
  0.6), not only against counts or ranges.
* Plot macros use `program cif` and `run`, never `program jmol` or
  `run plot:inter`.
* Add every test to `test_all.mac` and `manifest.toml` together. Tag it
  (`fast`, `mc`, …), and check with `tools/run_legacy_tests.py --check-test-all`.
* After adding tests, re-run `tools/command_coverage.py` and update the
  appendix below.

---

## Appendix – command coverage tables (generated)

Output of `tools/command_coverage.py` for the run above.

| Program | Commands reached | Sub-commands / options reached |
|---|---|---|
| suite | 3 / 7 | – |
| discus | 28 / 61 | 136 / 402 |
| kuplot | 34 / 82 | 9 / 69 |
| diffev | 0 / 28 | – |
| refine | 6 / 9 | 3 / 4 |
| experi | 0 / 1 | 0 / 10 |
| lib | 19 / 31 | – |

#### suite

**Not reached:** `diffev`, `parallel`, `show`, `update`

| Command | Tests | Sub-commands reached | Sub-commands not reached |
|---|---:|---|---|
| `discus` | 30 | – | – |
| `kuplot` | 6 | – | – |
| `refine` | 1 | – | – |

#### discus

**Not reached:** `anis`, `append`, `asym`, `bragg`, `change`, `copy`, `d2r`, `define`, `demolecularize`, `diff-four`, `displacement`, `exp2pdf`, `export`, `import`, `index2cell`, `inverse`, `kick`, `mc`, `molecularize`, `patterson`, `perioditize`, `private`, `project`, `r2d`, `remove`, `rmc`, `seed`, `spacegroup`, `storage`, `superspace`, `switch`, `ther`, `wyckoff`

| Command | Tests | Sub-commands reached | Sub-commands not reached |
|---|---:|---|---|
| `boundary` | 1 | – | – |
| `branch` | 25 | – | – |
| `chemistry` | 16 | 7 / 18 | angle, apply_symmetry, bval, coordination, disp, env, field, homo, mode, neig, trans |
| `connectivity` | 1 | 0 / 9 | add, create, delete, recreate, remove, reset, set, show, switch |
| `decorate` | 2 | 5 / 5 | – |
| `domain` | 26 | 7 / 8 | instructions |
| `find` | 2 | – | – |
| `fourier` | 13 | 21 / 27 | calc, delf, electron, energy, hkl, layer |
| `insert` | 1 | 0 / 14 | adp, cent, character, density, file, fuzzy, origin, reset, run, show, type, xaxis, yaxis, zaxis |
| `mmc` | 5 | 4 / 5 | apply_symmetry |
| `output` | 16 | 4 / 9 | input, range, reset, show, threshold |
| `pdf` | 1 | 5 / 9 | data, i/jdese, i/jsele, run |
| `plot` | 26 | 16 / 28 | asym, back, bonds, chem, col, hkl, line, mdes, mole, msel, title, type |
| `powder` | 3 | 5 / 7 | electron, neut |
| `property` | 3 | 1 / 2 | set/clear |
| `purge` | 6 | – | – |
| `read` | 30 | 4 / 5 | exit |
| `replace` | 11 | – | – |
| `reset` | 4 | – | – |
| `save` | 13 | 4 / 13 | asym, chem, des, format, incl, omit, prop, reset, show |
| `shear` | 2 | 6 / 20 | calc, des, eigen, hkl, mdes, mincl, mode, msel, odes, oincl, osel, show, uvw, vector |
| `show` | 11 | 5 / 18 | adp, asym, bval, chem, config, connectivity, envi, error, mass, menv, object, scat, symmetry |
| `stack` | 4 | 11 / 19 | ccol, cele, dist, four, random, reset, rotate, sigma |
| `surface` | 5 | 2 / 4 | character, show |
| `symmetry` | 3 | 14 / 24 | asym, calc, chem, des, hkl, mdes, odes, oincl, osel, use |
| `transform` | 1 | 0 / 33 | anew, aold, asnew, asold, asym, bnew, bold, bsnew, bsold, c2new, c2old, chem, cnew, cold, csnew, csold, des, h2new, h2old, incl, onew, oold, reset, run, sel, set, show, xnew, xold, ynew, yold, znew, zold |
| `vprod` | 1 | – | – |
| `waves` | 2 | 15 / 26 | asym, chem, dens, mdes, mrep, msel, opti, osci, rot, shift, tran |

#### kuplot

**Not reached:** `adt`, `alloc`, `bond`, `branch`, `close`, `color`, `convolute`, `cost`, `derivative`, `dsave`, `ecolor`, `etype`, `excl`, `fft`, `fill`, `fit`, `fourier`, `function`, `glat`, `grid`, `hcolor`, `hlab`, `hpak`, `htype`, `ident`, `integral`, `kpar_par`, `kpara`, `lart`, `mean`, `merge`, `mouse`, `nexus`, `print`, `ptype`, `rdef`, `rebin`, `sann`, `scale`, `sdef`, `show`, `sleg`, `sm`, `smax`, `sort`, `spline`, `torefine`, `window`

| Command | Tests | Sub-commands reached | Sub-commands not reached |
|---|---:|---|---|
| `achx` | 23 | – | – |
| `achy` | 23 | – | – |
| `achz` | 9 | – | – |
| `angle` | 10 | – | – |
| `aver` | 19 | – | – |
| `buffer` | 12 | – | – |
| `ccalc` | 1 | – | – |
| `cmap` | 9 | – | – |
| `fname` | 22 | – | – |
| `font` | 5 | 1 / 4 | color, just, typ |
| `frames` | 12 | 4 / 7 | bframe, cframe, frame |
| `fset` | 12 | – | – |
| `hart` | 9 | – | – |
| `hline` | 9 | – | – |
| `kcalc` | 1 | – | – |
| `ksave` | 1 | 1 / 6 | form, m999, run, show, threshold |
| `lcolor` | 3 | – | – |
| `load` | 28 | 3 / 16 | csv, de, gs, h5, ma, mc, mp, mrc, nexus, pg, sc/st, xx, zz |
| `ltype` | 2 | – | – |
| `lwidth` | 1 | – | – |
| `mark` | 24 | – | – |
| `match` | 2 | – | – |
| `mcolor` | 1 | – | – |
| `msize` | 19 | – | – |
| `mtype` | 2 | – | – |
| `orientation` | 3 | – | – |
| `plot` | 26 | – | – |
| `reset` | 26 | – | – |
| `rvalue` | 16 | – | – |
| `save` | 3 | – | – |
| `skal` | 24 | – | – |
| `smooth` | 1 | – | – |
| `tit1` | 22 | – | – |
| `tit2` | 22 | – | – |

#### diffev

**Not reached:** `adapt`, `allocate`, `backup`, `branch`, `compare`, `constraint`, `data`, `deallocate`, `dismiss`, `donor`, `fix`, `initialize`, `lastfile`, `logfile`, `newparam`, `purge`, `read`, `refine`, `release`, `reset`, `restart`, `restrial`, `run_mpi`, `selection`, `show`, `trialfile`, `type`, `write`

#### refine

**Not reached:** `finished`, `fix`, `sigma`

| Command | Tests | Sub-commands reached | Sub-commands not reached |
|---|---:|---|---|
| `data` | 1 | – | – |
| `newparam` | 1 | – | – |
| `reset` | 1 | – | – |
| `run` | 1 | – | – |
| `set` | 1 | 3 / 4 | log |
| `show` | 1 | – | – |

#### experi

**Not reached:** `exp2pdf`

#### lib

**Not reached:** `cd`, `fend`, `fexist`, `fsub`, `input`, `learn`, `lend`, `mount`, `seed`, `sleep`, `timer`, `umount`

| Command | Tests | Sub-commands reached | Sub-commands not reached |
|---|---:|---|---|
| `break` | 1 | – | – |
| `continue` | 10 | – | – |
| `do` | 24 | – | – |
| `echo` | 30 | – | – |
| `evaluate` | 16 | – | – |
| `exit` | 31 | – | – |
| `fclose` | 31 | – | – |
| `fformat` | 5 | – | – |
| `fget` | 1 | – | – |
| `fopen` | 31 | – | – |
| `fput` | 31 | – | – |
| `if` | 30 | – | – |
| `matrix` | 1 | – | – |
| `set` | 31 | – | – |
| `show` | 1 | – | – |
| `stop` | 30 | – | – |
| `system` | 30 | – | – |
| `variable` | 31 | – | – |
| `wait` | 19 | – | – |

#### `set` parameters used inside menus

| Menu | Parameters (number of tests) |
|---|---|
| discus/chemistry | bang (1), bin (2), blen (2), crystal (1), mode (7), nei (1), neig (1), vec (1) |
| discus/decorate | citric (1), water (2) |
| discus/domain | remove (1) |
| discus/fourier | aver (12), technique (1) |
| discus/mmc | allowed (1), cyc (3), cycl (2), feed (5), finish (4), mode (5), move (2), neig (5), output (1), targ (4), target (1), temp (5), vec (5) |
| discus/pdf | bound (1), corrlin (1), corrquad (1), density (1), finite (1), qbroad (1), qdamp (1), qmax (1), rad (1), rang (1), srat (1), therm (1), weight (1) |
| discus/plot | 0 (17), 1 (19), 2 (3), 3 (2), d001 (1), d010 (1), d100 (1) |
| discus/powder | axis (3), calc (3), delta (3), disp (3), dq (3), four (3), lpcor (3), profile (3), qmax (3), qmin (3), temp (3), wvle (3) |
| discus/stack | modulus (4), trans (4) |
| discus/surface | distance (3) |
