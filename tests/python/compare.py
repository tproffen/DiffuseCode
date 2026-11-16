"""Numeric comparisons and golden reference data for the Python tests.

Golden data live in ``tests/golden/<test id>/``: one ``<name>.npz`` per
result plus ``meta.json`` (git commit, date, platform, NumPy version).

To (re)create them, run the tests with ``DIFFUSE_UPDATE_GOLDEN=1``. Those
tests then write their result and are skipped instead of compared. Review
the diff before committing new golden files.
"""

from __future__ import annotations

import datetime as _dt
import json
import os
import platform
import subprocess
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[2]
GOLDEN = REPO / "tests" / "golden"


def rvalue(obs, calc) -> float:
    """Unweighted R-value sqrt(sum (o - c)^2 / sum o^2), as KUPLOT ``rval ..., one``."""
    obs = np.asarray(obs, dtype=np.float64)
    calc = np.asarray(calc, dtype=np.float64)
    if obs.shape != calc.shape:
        raise ValueError(f"shape mismatch: {obs.shape} vs {calc.shape}")
    return float(np.sqrt(np.sum((obs - calc) ** 2) / np.sum(obs ** 2)))


def load_columns(path: Path, ncol: int = 2) -> np.ndarray:
    """Read a whitespace-separated text file, skipping ``#`` header lines."""
    return np.loadtxt(path, comments="#", usecols=range(ncol), ndmin=2)


def _meta() -> dict:
    try:
        commit = subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"],
                                capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        commit = "unknown"
    return {
        "git_commit": commit,
        "date": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "numpy": np.__version__,
    }


def updating() -> bool:
    return os.environ.get("DIFFUSE_UPDATE_GOLDEN", "") not in ("", "0")


def check_golden(test_id: str, name: str, *, r_max: float, **arrays) -> None:
    """Compare ``arrays`` with ``golden/<test_id>/<name>.npz``.

    Axis arrays (any key except the last) must match to float32 rounding;
    the last array is the result and must have an R-value below ``r_max``.
    With DIFFUSE_UPDATE_GOLDEN=1 the file is (re)written and the test skipped.
    """
    folder = GOLDEN / test_id
    path = folder / f"{name}.npz"
    if updating():
        folder.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(path, **arrays)
        meta_path = folder / "meta.json"
        meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}
        meta[name] = _meta() | {"r_max": r_max}
        meta_path.write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n")
        pytest.skip(f"golden data written to {path.relative_to(REPO)}")
    if not path.exists():
        pytest.fail(f"no golden data {path.relative_to(REPO)}; "
                    "run once with DIFFUSE_UPDATE_GOLDEN=1")
    ref = np.load(path)
    assert sorted(ref.files) == sorted(arrays), f"golden keys {ref.files} != {list(arrays)}"
    *axes, result = arrays
    for key in axes:
        np.testing.assert_allclose(arrays[key], ref[key], rtol=1e-6, atol=1e-6,
                                   err_msg=f"axis '{key}' differs from golden")
    r = rvalue(ref[result], arrays[result])
    assert r < r_max, f"R-value of '{result}' against golden is {r:.3g} (limit {r_max})"
