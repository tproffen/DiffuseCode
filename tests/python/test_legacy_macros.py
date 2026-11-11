"""pytest view of the legacy macro tests in tests/manifest.toml.

The work is done by tools/run_legacy_tests.py; this file only makes each
manifest entry a pytest test, so the macro tests can sit next to future
Python-API tests.

    pytest tests/python                       # manifest checks only, if no binary
    DISCUS_SUITE=/path/to/discus_suite pytest tests/python -m "not slow"
    DISCUS_SUITE=... pytest tests/python -k "04."

discus_suite is found via $DISCUS_SUITE or PATH. Without it, the macro tests
are skipped and only the manifest checks run. Tests run one at a time here;
use the runner directly (-j N) for parallel runs.
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import sys
import tempfile
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
_spec = importlib.util.spec_from_file_location(
    "run_legacy_tests", REPO / "tools" / "run_legacy_tests.py")
runner = importlib.util.module_from_spec(_spec)
sys.modules["run_legacy_tests"] = runner      # needed by dataclasses
_spec.loader.exec_module(runner)

TESTS = runner.load_manifest(runner.DEFAULT_MANIFEST)
EXE = runner.find_exe(os.environ.get("DISCUS_SUITE"))


def test_manifest_matches_test_all():
    """Every macro that test_all.mac checks is in the manifest, and back."""
    assert runner.check_test_all(TESTS) == 0


def test_manifest_ids_unique_and_files_exist():
    # load_manifest raises on duplicates and on missing macro files
    assert len({t.id for t in TESTS}) == len(TESTS)


def _param(t):
    marks = [getattr(pytest.mark, tag) for tag in t.tags]
    if not t.enabled:
        marks.append(pytest.mark.skip(reason=t.skip_reason))
    if "mpi" in t.tags:
        marks.append(pytest.mark.skip(reason="MPI tests are not run from pytest"))
    return pytest.param(t, id=t.id, marks=marks)


@pytest.fixture(scope="session")
def workdir():
    path = os.environ.get("DISCUS_TEST_WORKDIR")
    return Path(path) if path else Path(tempfile.mkdtemp(prefix="discus_pytest_"))


@pytest.mark.skipif(EXE is None, reason="discus_suite not found ($DISCUS_SUITE or PATH)")
@pytest.mark.parametrize("test", [_param(t) for t in TESTS])
def test_macro(test, workdir):
    opts = argparse.Namespace(
        xvfb=not os.environ.get("DISPLAY") and runner.shutil.which("xvfb-run") is not None,
        omp_threads=None, timeout_scale=1.0, keep=False)
    displays = runner.queue.Queue()
    displays.put(100 + (os.getpid() % 100) * 8)
    r = runner.run_one(test, EXE, workdir, displays, opts)
    if r.status == runner.RAN:
        pytest.skip(f"ran, but has no check: {r.message}")
    assert r.status == runner.PASS, f"{r.status}: {r.message}\nlog: {r.log}"
