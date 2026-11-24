"""Register the manifest tags as pytest markers (tests/manifest.toml)."""

TAGS = {
    "fast": "runs in under 10 s",
    "medium": "runs in under 30 s",
    "slow": "runs for 30 s or more",
    "mc": "stochastic (Monte Carlo / random); results vary between runs",
    "mpi": "needs MPI",
    "nocheck": "has no self-check yet",
    "exclusive": "sets its own OpenMP thread count; do not run in parallel",
}


def pytest_configure(config):
    for tag, text in TAGS.items():
        config.addinivalue_line("markers", f"{tag}: {text}")


def pytest_report_header(config):
    """Say which discus_suite the legacy macro tests use (test_legacy_macros.py)."""
    import importlib.util
    from pathlib import Path
    tool = Path(__file__).resolve().parents[2] / "tools" / "run_legacy_tests.py"
    spec = importlib.util.spec_from_file_location("_legacy_runner_header", tool)
    runner = importlib.util.module_from_spec(spec)
    import sys
    sys.modules[spec.name] = runner              # needed by dataclasses
    spec.loader.exec_module(runner)
    exe, where = runner.find_exe_where(None)
    if exe is None:
        return "discus_suite: not found (legacy macro tests skipped; set DISCUS_SUITE)"
    return runner.describe_exe(exe, where)


# ------------------------------------------------------- diffuse
import os
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))                    # compare.py

# Which diffuse package to test: with DIFFUSE_PYTHON_PATH (a CMake build
# folder) the source in python/; otherwise an installed one ("pip install
# ."), and the source only if none is installed.
import importlib.util  # noqa: E402

SOURCE_PACKAGE = HERE.parents[1] / "python"
USE_SOURCE = bool(os.environ.get("DIFFUSE_PYTHON_PATH")) or \
    importlib.util.find_spec("diffuse") is None
if USE_SOURCE:
    sys.path.insert(0, str(SOURCE_PACKAGE))


@pytest.fixture(scope="session")
def suite():
    """The in-process DISCUS suite (``diffuse.session()``).

    Skips the test if the compiled core (libdiffuse_core) cannot be
    loaded. Build
    it with -DDIFFUSE_PYTHON=ON and set DIFFUSE_PYTHON_PATH=<build>/python
    (DIFFUSE_REQUIRE_PYTHON=1 turns the skip into an error, as in CI).
    The Fortran state is global, so there is one session per pytest run.
    """
    import diffuse
    try:
        return diffuse.session()
    except ImportError as exc:
        if os.environ.get("DIFFUSE_REQUIRE_PYTHON"):
            raise
        pytest.skip(str(exc).splitlines()[0])
