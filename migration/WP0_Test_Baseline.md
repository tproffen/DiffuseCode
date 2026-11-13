# WP0 – Regression test baseline

Part of the Python migration plan ([`DiffuseCode_Python_Migration_Plan.md`](DiffuseCode_Python_Migration_Plan.md), option B).
WP0 is the safety net: before any Fortran is wrapped or refactored, we need a test
suite that runs unattended, reports pass/fail reliably, and stores reference
("golden") results to compare against.

This document has three parts:

1. **Where we are** – what the existing `tests/` folder covers and what happened
   when it was run headless on Linux.
2. **WP0 work plan** – the tasks to turn it into an automated regression baseline.
3. **Progress log** – what has been done so far.

`mixscat` is deprecated and is left out of WP0 (and of the migration).

---

## 1. Where we are

### 1.1 What is in `tests/`

| Folder | Topic | `.mac` files | Macros with a self-check | Reference data (`RESULTS/`) | Run by `test_all.mac` |
|---|---|---:|---:|---|---|
| 00.Introduction | Command-language basics | 5 | 0 | – | no |
| 01.Simulation | Building structures, symmetry, shear | 19 | 6 | – | 4 tests |
| 02.Calculation_Scattering | 1D/2D Fourier, lots, powder | 29 | 5 | yes | 5 tests |
| 03.Pair_Distribution_Function | PDF | 11 | **0** | – | runs, **no check** |
| 04.SRO | Short-range order via MMC | 138 | 6 | yes | 4 tests (SRO_2 commented out) |
| 05.Stack | Stacking faults | 48 | 4 | yes | 4 tests |
| 06.Domains | Domains | 32 | 4 | yes | 4 tests |
| 07.Refinement | DIFFEV/REFINE | 50 | 1 | – | 1 test |
| 08.Nanoparticle_example | Nanoparticles, powder, PDF | 41 | 3 | yes | 2 tests |
| 09.Nanoparticle_refinement | MPI refinement | 364 | 0 | – | 1 shell script, **does not run** (see 1.3) |
| 11.Properties | Atom properties | 7 | 1 | – | 1 test |
| 12.Waves | Modulation waves | 8 | 3 | yes | 2 tests |
| 13.Analyze | Neighbour/bond analysis | 10 | 1 | yes | 1 test |
| 14.Decoration | Surface decoration | 46 | 10 | – | 2 tests |
| 15.Timing | OpenMP benchmarks | 15 | 0 | – | no |

**How a test works.** `test_all.mac` changes into each folder, runs a macro with
parameters, then calls `check.mac`. Each macro sets `LTEST_MAC = TRUE` if its own
check passes. A check is either:

* a geometric value (`blen`, `bang`, atom counts), or
* a KUPLOT comparison against a stored file: `load` the new result and
  `RESULTS/<file>`, then `rval` (R-value) below a threshold, e.g. `< 0.001` for
  Fourier and powder, `< 0.25` for the Monte-Carlo SRO examples.

Total: **30 test invocations**, about 3.5 minutes on 4 cores.

### 1.2 Headless run results (2026-11-04)

Build: Ubuntu 24.04, gfortran 13.3, CMake 3.28, FFTW, HDF5, `pgplot5` from the
Ubuntu repository, run under `xvfb-run` (virtual X display), no human input.

| Run | Result |
|---|---|
| `test_all.mac` as shipped | Stopped in 05.Stack after 4 sections: `***APPL*** Interactive plot for JMOL only (-152)`. **Process exit code was 0.** |
| Each test run in isolation | 22 PASS, 7 ABORTED (all the same JMOL error), 1 ran without a check |
| With JMOL lines disabled (this branch) | **All 13 sections SUCCESS, 200 s, exit code 0.** CIF output files are still written; `SIMPLE/PLOT/block.cif` is byte-identical to the committed copy. |

**JMOL fix in this branch.** `program jmol` first runs `which jmol`. If JMOL is
not installed it raises error -152, which aborts the whole macro run.
`run plot:inter` also raises -152 for anything other than JMOL. In the 9 plot
macros that `test_all.mac` reaches, both lines are now commented out and replaced
with `program cif` + `run`. These write the same CIF file through the same code
path (`discus_plot.f90:1448`), just without opening the viewer.

Files changed:
`05.Stack/STACK/plot.mac`, `05.Stack/STACK_FILE/plot.mac`,
`06.Domains/CUBE/plot.cif.mac`, `08.Nanoparticle_example/{SIMPLE,ELLIPSOID,Q}/plot.cif.mac`,
`14.Decoration/Deco_01/plot.mac`, `14.Decoration/Deco_02/plot.mac`,
`14.Decoration/Deco_02/clean_surface.mac`.

Other tutorial macros still call JMOL interactively, but the suite does not run
them: `01.Simulation/plotcif.mac`, `04.SRO/SRO_MOLECULES/plot_jmol.mac`,
`05.Stack/{STACK_DISTRIBUTION,WURTZITE}`, `06.Domains/Molecules`,
`08.Nanoparticle_example/COMPLEX/main.mac`, `11.Properties/{prop_ex_2,plot.cif}.mac`,
`15.Timing/OpenMP/plot.mac`.

### 1.3 Problems found

**Test-harness problems**

1. **Errors are invisible to automation.** A macro error ends `discus_suite` with
   exit code 0, and so does a failing check (`stop` in `check.mac`). The only
   pass/fail signal is text on screen ("SUCCESS in …" / "Failure, macro …").
   With macro echo on, the log also contains the *source* of both branches of
   every `if`, so grepping for "Failure" finds false hits.
2. **One failure hides all later tests.** `test_all.mac` is a single sequential
   macro, and the first error ends the run.
3. **Some sections report SUCCESS without testing anything:**
   * 03.Pair_Distribution_Function – `pdf.mac` runs but has no check.
   * 09.Nanoparticle_refinement – `system ./TEST.sh` fails with
     *Permission denied* (no execute bit in git), and needs `mpiexec` and
     `DISCUS_NCPU`. The section still prints SUCCESS.
   * 14.Decoration/Deco_01 – `all.mac` has no aggregated check.
4. **Randomness is not controlled.** The suite has a `seed` command
   (`seed <value>` or `seed <a>,<b>,…, group:<n>`, see `lib_f90/lib_f90.hlp`).
   Only 12 macros use it (mostly `diffev_best.mac`), and none of the tests in
   `test_all.mac` do. The checks therefore use loose tolerances (SRO
   correlation ±0.08, R < 0.25).
   **`seed` does not currently make MMC runs reproducible.** A quick check ran
   `04.SRO/SRO_0/sro.0.mac 0.7,1.5` twice with `seed 1234`, both at the
   suite prompt and directly before `mmc`, with `OMP_NUM_THREADS=1` and 4.
   Every time the two runs gave different MC acceptance counts, and about
   half the lines of `example0.stru` differed. This needs debugging (task
   0.5a) before golden MC data can be exact.
5. **Generated output is committed next to inputs.** Example:
   `02.Calculation_Scattering/example3.inte` sits beside
   `RESULTS/example3.inte`. `tests/` is 109 MB / 2,053 files, and it is hard
   to tell inputs, references and leftovers apart.
6. **Coverage gaps.** A rough grep finds no use of these DISCUS commands in
   the test macros: `anis, append, asym, bragg, change, define, copy,
   demolecularize, interpret, inverse, kick, molecularize, patterson,
   perioditize, rmc, spacegroup, storage, super, switch, wyckoff`. Some may be
   reached indirectly or under abbreviations; this needs confirming.

**Build problems**

7. `DIFFUSE_GRAPHICS=OFF` does not work: `suite/prog/CMakeLists.txt:5`
   always calls `find_package(PGPLOT REQUIRED)`.
8. `make install` assumes PGPLOT's `grfont.dat` and `pgxwin_server` are in
   `PGPLOT_DIR`. On Debian/Ubuntu they are in `/usr/lib/pgplot5` and
   `/usr/bin`. The binary works fine from the build tree.
9. The build requires an `FC` environment variable even for re-running CMake.

**Found while doing 0.2 / 0.3**

10. **Start-up race (fixed).** On every start `discus_suite` pings GitHub
    (`check_github` in `lib_f90/appl_unix.f90`) and wrote the answer to a
    fixed file, `/tmp/check_github.txt`. Two processes starting at the same
    time could delete each other's file, and one crashed with a Fortran
    runtime error. The file name now includes the process ID, like the other
    temporary files. The ping also means every start tries to reach the
    network; in CI it simply fails quickly.
11. **Macros override the thread count.** `04.SRO/SRO_1/sro.1.mac` and
    `SRO_3/sro.3.mac` call `set parallel, useomp:use, nthreads:all`, which
    overrides `OMP_NUM_THREADS`. Next to other tests this over-commits the
    CPUs, and the OpenMP threads slow down by up to ~100× (2.4 s → 130–340 s).
    These tests are tagged `exclusive` and the runner runs them alone.
12. **Endless prompt loop.** When `check.mac` calls `stop` (interactive mode)
    and stdin is at end of file, `discus_suite` prints `suite >` prompts
    forever (18 MB of log in 20 s). Batch mode never calls `stop`, and the
    runner's timeout would catch it, but the input loop should end on EOF.

---

## 2. WP0 work plan

Goal: one command (and one CI job) that builds DiffuseCode, runs every test in
isolation, writes machine-readable results, and **fails** when a test fails.
It also stores golden numerical outputs that later Python-API tests (WP3–WP8)
can compare against.

### Proposed layout

```
tests/                       # existing macro tests stay where they are
  check.mac                  # [done] batch mode writes PASS/FAIL to a result file
  test_all.mac               # kept for interactive use (unchanged)
  manifest.toml              # [done] list of tests (TOML: read by Python's standard library)
  README.md                  # [done] how to run the tests
  golden/                    # 0.6: reference outputs, one folder per test id
  python/                    # pytest harness
    conftest.py              #   [done] registers the manifest tags as markers
    test_legacy_macros.py    #   [done] one test per manifest entry + manifest checks
    compare.py               #   0.6: numeric comparisons (R-value, max-abs, atom sets)
tools/
  run_legacy_tests.py        # [done] runner -> console table, JUnit XML, JSON
.github/workflows/
  legacy-tests.yml           # 0.9: build + run on every push to claude-python-app
```

### Tasks

| # | Task | Details | Done when | Effort |
|---|---|---|---|---|
| 0.1 | **Headless run** | JMOL disabled in the test plot macros *(done in this branch)*. Run under `xvfb-run` with PGPLOT; later use the no-graphics build once 0.8 is fixed. | `test_all.mac` completes unattended | done |
| 0.2 ✅ | **Machine-readable results** | Extend `check.mac` to append `<dir> <macro> PASS/FAIL` to `test_results.txt` instead of `stop`-ing. Add a driver that runs each manifest entry in a fresh temp copy (a prototype exists from this analysis) and converts results to JUnit XML. Treat a missing result line as ABORTED. | CI shows per-test pass/fail; process exits non-zero on any failure or abort | 2–3 d |
| 0.3 ✅ | **Test manifest** | Move the list of tests out of `test_all.mac` into `manifest.toml`: test id, folder, macro, args, expected outputs, comparison type and tolerance, and tags (`fast`, `mc`, `mpi`, `slow`). | Harness and `test_all.mac` cover the same 30 tests | 1–2 d |
| 0.4 | **Close the silent gaps** | Add a check to 03 PDF (compare G(r) against a new `RESULTS/` file). Re-enable or remove SRO_2. Add an aggregated check to Deco_01. Make 09 runnable: execute bit + serial variant for CI, MPI variant tagged `mpi`. | No section prints SUCCESS without a check | 2–3 d |
| 0.5a | **Debug `seed`** | Find out why `seed` does not give repeatable MMC results (see 1.3 #4). Candidates: the generator being re-seeded inside `mmc`/`read`/`insert`, per-thread generators, or a seed state that is not fully applied. Add a minimal test: seed → draw N numbers → seed again → same N numbers, then the same through `mmc run`. | Two runs with the same seed give identical `sro.0` output on one build | 1–3 d |
| 0.5 | **Reproducibility** | Once 0.5a is fixed, add `seed` to every stochastic test (MC, RMC, DIFFEV, stacking, domains). Record whether a test is bit-reproducible on one platform/thread count, and keep statistical tolerances where it is not (OpenMP order, compiler). | Same seed gives the same output on the same build | 2 d |
| 0.6 | **Golden data capture** | For each test, copy key outputs (intensities, powder, PDF, 3D-PDF, final structure/CIF, histograms) into `tests/golden/<id>/` with a `meta.json` (DiffuseCode version, git commit, compiler, flags, date). Python compares with NumPy (R-value, max-abs, atom-set match) using tolerances from the manifest. Keep the macro self-checks as a second line of defence. | Every test has at least one golden numeric output | 3–4 d |
| 0.7 | **Coverage expansion** | Confirm the uncovered-command list from 1.3, then add short, deterministic tests for the commands that WP3–WP8 will wrap first: space-group expansion / `wyckoff`, `patterson`, `inverse`, `rmc`, `molecularize`/`demolecularize`, `perioditize`, `storage`/`super`, `bragg`, and kuplot `fit`. | All commands planned for WP3–WP8 have at least one test | 5–8 d |
| 0.8 | **Build fixes** | Make `suite/prog/CMakeLists.txt` honour `DIFFUSE_GRAPHICS`. Make the PGPLOT install paths configurable / optional. Drop the hard `FC` requirement in favour of CMake's normal compiler detection. | `cmake -DDIFFUSE_GRAPHICS=OFF` builds and the suite passes | 1–2 d |
| 0.9 | **CI** | GitHub Actions on `ubuntu-24.04`: apt (gfortran, FFTW, HDF5, readline, PNG, X11, pgplot5, xvfb), build, run `tools/run_legacy_tests.py --tags fast,mc`, upload JUnit + logs. Nightly job for `slow`/`mpi`. Add a macOS job once Linux is stable. | Red/green status on every push to `claude-python-app` | 1–2 d |
| 0.10 | **Repo hygiene** | Run `cleanup.mac` logic, remove committed generated outputs, add a `.gitignore` for them, and keep only inputs + `RESULTS/` + `golden/`. Consider Git LFS for golden files larger than ~1 MB. | `git status` is clean after a test run | 1 d |

**Total: about 4–6 weeks for one developer.** That fits the "months 0–3"
slot in the plan, with room for coverage work alongside WP1.

### Suggested order

1. **0.2 + 0.3** – reliable pass/fail first. Everything else depends on it.
2. **0.9** – CI right away, so regressions are seen from day one.
3. **0.4 + 0.5a + 0.5** – remove false SUCCESS, fix `seed`, and make results reproducible.
4. **0.6** – golden data, the input for the Python-API work.
5. **0.8 + 0.10** – build and repo clean-up, then turn CI to the no-graphics build.
6. **0.7** – coverage growth, continuing in parallel with WP1–WP3.

### Open questions for the team

* Is changing `check.mac`'s behaviour (write a result instead of `stop`) OK for
  interactive users of `test_all.mac`? An `LTEST_BATCH` variable could switch
  between the two.
* Which tolerances count as "same physics" for MC-based tests, and who signs
  off on them?
* Should golden data live in the repo (Git LFS) or in a separate data repository
  or release asset?
* Is MPI coverage needed in CI, or is a nightly or manual run enough?

### Reproducing the headless run

```bash
sudo apt-get install -y gfortran cmake libfftw3-dev libhdf5-dev libreadline-dev \
    libpng-dev libx11-dev zlib1g-dev pgplot5 xvfb
```

```bash
mkdir -p ../DiffuseBuild && cd ../DiffuseBuild && \
FC=$(which gfortran) cmake ../DiffuseCode -DCMAKE_BUILD_TYPE=Release -DPGPLOT_DIR=/usr/lib/pgplot5 && \
make -j"$(nproc)"
```

Unattended, with per-test results (from the `DiffuseCode` folder; the runner
copies each test, so `tests/` stays clean):

```bash
python3 tools/run_legacy_tests.py --exe ../DiffuseBuild/suite/prog/discus_suite -j 4 --junit results.xml
```

Or the classic way, on a copy because the macros write output files into
`tests/`:

```bash
cp -r tests /tmp/dc_tests && cd /tmp/dc_tests && \
PGPLOT_DIR=/usr/lib/pgplot5 xvfb-run -a "$OLDPWD/../DiffuseBuild/suite/prog/discus_suite" -macro test_all.mac < /dev/null
```

---

## 3. Progress log

### 0.2 + 0.3 – runner, batch mode and manifest

* `tests/check.mac` has a batch mode. With `LTEST_BATCH = TRUE` it appends
  `PASS <name>` / `FAIL <name>` to the file in `LTEST_FILE` and carries on,
  judging each check on its own. Interactive behaviour is unchanged
  (`test_all.mac` still passes all 13 sections). `check.mac` declares
  `LTEST_BATCH` itself, so test macros run by hand still work.
* `tests/manifest.toml` lists 33 tests: the 31 macros run by `test_all.mac`,
  plus `04.sro2` and `09.znse_pdf_mpi`, which are disabled with a reason.
  `--check-test-all` confirms that it matches `test_all.mac`.
  `14.deco01_all` now counts its 7 internal checks; before, it had no
  overall check.
* `tools/run_legacy_tests.py` (standard library only):
  * Runs each test in a fresh copy, in its own process group, under xvfb when
    there is no display. Runs tests in parallel (`-j`) and runs tests tagged
    `exclusive` alone.
  * Kills the whole process group on timeout, retries on another display if
    Xvfb cannot start, and writes JUnit XML and JSON.
  * Exits with 1 on FAIL, ABORTED or TIMEOUT.
* `tests/python/` gives the same tests as pytest tests, plus the manifest
  checks (which run without a binary).
* `tests/README.md` explains how to run and add tests.
* `lib_f90/appl_unix.f90`: fixed the start-up race (finding 10).

Checked:

| Run | Result |
|---|---|
| Full manifest, `-j 4` (4 cores) | 30 PASS, 1 RAN (`03.pdf`), 2 SKIPPED; 78 s wall time; repeated 3 times, same result |
| Full manifest, `-j 1` | same result; 210 s |
| Deliberate failures (wrong result, bad parameter, too few checks, 5 s timeout) | FAIL, ABORTED (first error line reported), ABORTED, TIMEOUT; exit code 1; no processes left behind |
| Display already taken | retried on another display, PASS |
| Runner and pytest at the same time | all pass; no display clashes |
| `test_all.mac` (interactive mode) | all 13 sections SUCCESS |
