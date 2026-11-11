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
| `--list` | show what would run |
| `--check-test-all` | check that `manifest.toml` and `test_all.mac` list the same macros |

Without an X display it runs each test under `xvfb-run` (KUPLOT needs one).
On Debian/Ubuntu install `xvfb` and `pgplot5`.

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

## Adding a test

1. Write the macro so that it sets `LTEST_MAC = TRUE` when its result is
   correct (see `02.Calculation_Scattering/1d.mac` for a comparison against
   `RESULTS/`).
2. Add it to `test_all.mac` (`@macro args` followed by `@../check.mac macro`).
3. Add a `[[test]]` entry to `manifest.toml`.
4. Run `python3 tools/run_legacy_tests.py --check-test-all` and the test itself
   with `--id`.

Plot macros must not need an interactive viewer: use `program cif` and a
plain `run` instead of `program jmol` and `run plot:inter`.
