# DISCUS-suite macro tests

The folders `00.Introduction` … `15.Timing` hold example and test macros.
Most of them check their own result and set `LTEST_MAC = TRUE` when it is
correct. `check.mac` evaluates that flag.

There are two ways to run the tests.

## 1. Interactively: `test_all.mac`

```
discus_suite -macro test_all.mac
```

Runs all checked macros in one session. The first failure prints
`Failure, macro … did not succeed` and stops. This is unchanged.

## 2. Unattended: `tools/run_legacy_tests.py`

For CI and for checking a build. It runs every test listed in
`manifest.toml` in its own fresh copy of the test folder and its own
`discus_suite` process. It reports PASS / FAIL / ABORTED / TIMEOUT per test
and exits with code 1 if any test did not pass.

```bash
python3 tools/run_legacy_tests.py --exe /path/to/discus_suite -j 4
```

Useful options:

| Option | Meaning |
|---|---|
| `--exe PATH` | the `discus_suite` binary (or set `$DISCUS_SUITE`, or have it on `PATH`) |
| `-j N` | run N tests in parallel (tests tagged `exclusive` still run alone, at the end) |
| `--tags fast` | only tests with one of these tags (`fast`, `medium`, `slow`, `mc`, …) |
| `--exclude-tags slow,mpi` | skip tests with these tags (default: `mpi`) |
| `--id '04.*'` | only tests whose id matches (repeatable) |
| `--junit FILE`, `--json FILE` | machine-readable reports |
| `--workdir DIR`, `--keep` | where copies and logs go; keep passing tests' copies too |
| `--seed N`, `--no-seed` | tests tagged `mc` run with `seed N` (default 1234; a test's own `seed` in the manifest wins) on one OpenMP thread, so they repeat exactly; `--no-seed` runs them unseeded with their own thread settings |
| `--retry-mc N` | re-run a failed test tagged `mc` up to N times; a pass on retry is reported as `flaky` (no longer used by CI: seeded `mc` tests repeat exactly) |
| `--list` | show what would run |
| `--check-test-all` | check that `manifest.toml` and `test_all.mac` list the same macros |

With a graphics build and no X display it runs each test under `xvfb-run`
(KUPLOT needs one); on Debian/Ubuntu install `xvfb` and `pgplot5`. A build
with `-DDIFFUSE_GRAPHICS=OFF` needs neither: use `--xvfb off`.

Results of a run with `-j 4` on 4 cores: 30 PASS, 1 RAN (`03.pdf`, no check
yet), 2 SKIPPED (`04.sro2`, `09.znse_pdf_mpi`), about 80 s wall time.

### How batch mode works

The runner writes a small driver macro that sets `LTEST_BATCH = TRUE` and
`LTEST_FILE` before calling the test macro. In batch mode `check.mac` appends
`PASS <name>` or `FAIL <name>` to `LTEST_FILE` instead of stopping. The driver
writes `END` as its last line.

| Result file | Status |
|---|---|
| expected number of `PASS` lines, then `END` | PASS |
| any `FAIL` line | FAIL |
| no `END`, or too few `PASS` lines | ABORTED (a command error ended the macro run; the first `***…***` error line is reported) |
| process killed after `timeout` | TIMEOUT |
| `check = "none"` and `END` | RAN (the macro finished, nothing was verified) |

`check.mac` declares `LTEST_BATCH` itself, so a test macro run by hand still
behaves as before.

### pytest

`tests/python/test_legacy_macros.py` turns each manifest entry into a pytest
test, and also checks the manifest against `test_all.mac`:

```bash
DISCUS_SUITE=/path/to/discus_suite pytest tests/python -m "not slow"
```

Without `DISCUS_SUITE` only the manifest checks run.

## 3. Python: `tests/python/test_suite_python.py`, `test_diffuse_api.py`

Run DISCUS in-process through the f2py module `suite_python` and check the
results as NumPy arrays, against the legacy `RESULTS/` files and against
golden data in `tests/golden/<id>/`. `test_diffuse_api.py` tests the
`diffuse` package (`python/diffuse/`, see `python/README.md`) and the demo
script `python/examples/powder_demo.py` (needs matplotlib). Build with `-DDIFFUSE_PYTHON=ON`, then:

```bash
pip install -r tests/requirements.txt      # gemmi, for independent references
DIFFUSE_PYTHON_PATH=/path/to/build/python pytest tests/python
```

Besides the API tests, `test_wp3_structures.py` and `test_wp4_fourier.py`
are the safety net for the WP3 (structures, I/O) and WP4 (single-crystal
Fourier) work: they compare DISCUS with gemmi's space groups and form
factors and with physical rules, and hold golden data. Tests marked
`xfail(strict=True)` document known DISCUS bugs; they turn red when the bug
is fixed, so the mark can be removed. Without gemmi the cross-checks are
skipped.

| Variable | Meaning |
|---|---|
| `DIFFUSE_PYTHON_PATH` | folder that holds `libdiffuse_core.so` and `suite_python*.so` (`<build>/python`) |
| `DIFFUSE_BACKEND` | `ctypes` (the C library `libdiffuse_core`, the default when it is there) or `f2py` (`suite_python`); CI runs the tests on both |
| `DIFFUSE_REQUIRE_PYTHON=1` | fail instead of skip if the module cannot be imported (used in CI) |
| `DIFFUSE_UPDATE_GOLDEN=1` | write the golden files instead of comparing; review the diff before committing |

To add a test, use the `suite` fixture (the session from `diffuse.session()`) or
the `diffuse` functions directly. `suite.run(commands, folder)` runs
macro text with `folder` as working directory and raises `MacroError` if a
command fails; `suite.value("n[1]")` evaluates an expression;
`suite.dataset(1)` / `suite.dataset_2d(1)` return KUPLOT data set 1. Load
results into KUPLOT (`kuplot`, `rese`, `load ...`) to read them back.
`compare.check_golden(id, name, r_max=..., axis..., result)` handles golden
data. The Fortran state is global and shared by all tests in one pytest run,
so start each calculation with the relevant `reset`.

## CI

`.github/workflows/tests.yml` builds with the Python module and runs the
pytest tests, then the legacy tests tagged `fast` and `mc`.

It runs only what the changes need (`tools/ci_changes.py`, first job
"changes"). A push to a `claude/**` branch is compared with
`claude-python-app`, a push to `claude-python-app` with the commit before:

| Changed | Runs |
|---|---|
| `migration/`, `*.md` outside `docs/`, `License.txt` | nothing |
| `docs/` | docs (Ubuntu 24.04) |
| Python package, `tests/python/`, `tests/golden/` | pytest (graphics build, Ubuntu 24.04 and 22.04), docs |
| Fortran sources, legacy tests, anything else | build, legacy tests, pytest (graphics build, 22.04) |
| `pyproject.toml`, `environment.yml`, CMake, `python/suite_python.f90`, `.github/` | everything, conda included |

Opening a pull request into `claude-python-app` runs nothing when the branch
already contains the latest `claude-python-app`: the push runs tested exactly
that code. If `claude-python-app` has moved on since the branch began, the
pull request run tests the merge result, for the jobs the PR's changes need;
that catches clashes with what was merged in between. Later pushes to the
branch run as pushes. A manual run also runs everything. GitHub shows the "Run workflow" button (Actions tab)
only for workflows on the default branch (`master`), so until the first
release start it from the command line, for any branch:
`gh workflow run tests.yml --ref claude-python-app`, or ask Claude. To see why a job ran or
not, open the "changes" job, or try it locally:
`python3 tools/ci_changes.py --files discus/prog/fourier.f90`.

## Notebooks in the repository

Notebooks (`python/examples/*.ipynb`) are committed without outputs; a test
in `tests/python/test_diffuse_api.py` checks it. To work on one in place,
set up [nbstripout](https://github.com/kynan/nbstripout) once per clone:

```bash
pip install nbstripout
nbstripout --install
git config filter.nbstripout.extrakeys "metadata.language_info"
```

Git then strips outputs and the run's Python details when comparing and
committing, so `git diff` stays empty after a run and results are never
committed. `git status` and `git pull` may still list the notebook (its size
changed); `git add` on it stores nothing new and clears that, and a pull then
replaces the file (the run's outputs are gone). To just run the examples,
use a copy of the folder instead.

## Adding a test

1. Write the macro so that it sets `LTEST_MAC = TRUE` when its result is
   correct (see `02.Calculation_Scattering/1d.mac` for a comparison against
   `RESULTS/`).
2. Add it to `test_all.mac` (`@macro args` followed by `@../check.mac macro`).
3. Add a `[[test]]` entry to `manifest.toml`.
4. Run `python3 tools/run_legacy_tests.py --check-test-all` and the test itself
   with `--id`.
5. Keep generated files out of git: references go into `RESULTS/` (or
   `tests/golden/` for the Python tests); everything the macro writes
   besides that belongs in `tests/.gitignore`. After running the test,
   `git status` should show only the files you meant to add.

Plot macros must not need an interactive viewer: use `program cif` and a
plain `run` instead of `program jmol` and `run plot:inter`.
