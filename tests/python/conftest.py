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


# ---------------------------------------------------------------- suite_python
import os
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))   # suite_session, compare


@pytest.fixture(scope="session")
def suite():
    """The in-process DISCUS suite (f2py module suite_python).

    Skips the test if the module cannot be imported. Build it with
    -DDIFFUSE_PYTHON=ON and set DIFFUSE_PYTHON_PATH=<build>/python.
    The Fortran state is global, so there is one session per pytest run.
    """
    from suite_session import SuiteSession, import_suite_python
    try:
        module = import_suite_python()
    except ImportError as exc:
        if os.environ.get("DIFFUSE_REQUIRE_PYTHON"):
            raise
        pytest.skip(f"suite_python not importable ({exc}); "
                    "set DIFFUSE_PYTHON_PATH=<build>/python")
    return SuiteSession(module)
