"""Tests for the demo Python API: diffuse.read_structure, diffuse.powder, plotting.

    DIFFUSE_PYTHON_PATH=<build>/python pytest tests/python/test_diffuse_api.py
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from compare import check_golden

REPO = Path(__file__).resolve().parents[2]
EXAMPLES = REPO / "python" / "examples"
ANATASE = EXAMPLES / "anatase.cif"
SCATTERING = REPO / "tests" / "02.Calculation_Scattering"


@pytest.fixture
def dc(suite):
    """The diffuse package, with the DISCUS session up."""
    import diffuse
    return diffuse


# ------------------------------------------------------------- structures
def test_read_cif_unit_cell(dc):
    cr = dc.read_structure(ANATASE)
    assert cr.natoms == 12                      # Ti4 O8 in the I41/amd cell
    assert cr.lattice == pytest.approx((3.73, 3.73, 9.37, 90, 90, 90))
    assert cr.volume == pytest.approx(3.73 * 3.73 * 9.37, rel=1e-6)
    names, xyz = cr.atoms()
    assert sorted(set(names)) == ["O2-", "TI4+"]
    assert (names == "TI4+").sum() == 4
    assert xyz.shape == (12, 3)
    assert np.all((xyz >= -0.5) & (xyz <= 1.0))


@pytest.mark.parametrize("supercell, natoms", [((2, 3, 1), 72), (2, 96)])
def test_read_cif_supercell(dc, supercell, natoms):
    cr = dc.read_structure(ANATASE, supercell=supercell)
    assert cr.natoms == natoms
    assert len(cr.positions()) == natoms


def test_read_stru_as_is(dc):
    cr = dc.read_structure(SCATTERING / "example3.stru")
    assert cr.natoms == 2500
    assert cr.supercell is None


def test_read_errors(dc, tmp_path):
    with pytest.raises(FileNotFoundError):
        dc.read_structure(tmp_path / "missing.cif")
    bad = tmp_path / "x.xyz"
    bad.write_text("1\n\nC 0 0 0\n")
    with pytest.raises(dc.DiffuseError):
        dc.read_structure(bad)
    with pytest.raises(ValueError):
        dc.read_structure(ANATASE, supercell=(0, 1, 1))


def test_working_directory_kept(dc, tmp_path, monkeypatch):
    """DISCUS's cd must not change the caller's working directory."""
    monkeypatch.chdir(tmp_path)
    dc.read_structure(ANATASE)
    assert Path.cwd() == tmp_path
    with pytest.raises(dc.MacroError):
        dc.run_macro("discus\n  no_such_command\nexit")
    assert Path.cwd() == tmp_path


def test_macro_error_message(dc):
    with pytest.raises(dc.MacroError) as info:
        dc.run_macro("discus\n  no_such_command 1\nexit")
    assert "Unknown command" in str(info.value)
    assert "***" in info.value.output


# ----------------------------------------------------------------- powder
def anatase_pattern(dc, **kw):
    cr = dc.read_structure(ANATASE, supercell=(10, 10, 4))
    return dc.powder(cr, qmin=0.5, qmax=6.0, dq=0.005, **kw)


def test_powder_grid_and_peaks(dc):
    p = anatase_pattern(dc)
    assert p.q[0] == pytest.approx(0.5) and p.q[-1] == pytest.approx(6.0)
    assert np.allclose(np.diff(p.q), 0.005)
    assert np.all(p.intensity >= 0)
    # strongest X-ray line of anatase is (101): Q = 2 pi sqrt(1/a^2 + 1/c^2)
    a, c = 3.73, 9.37
    q101 = 2 * np.pi * np.sqrt(1 / a**2 + 1 / c**2)
    assert p.q[np.argmax(p.intensity)] == pytest.approx(q101, abs=0.01)
    assert p.two_theta()[0] == pytest.approx(np.degrees(2 * np.arcsin(0.5 * 0.2 / (4 * np.pi))))
    check_golden("api.powder_anatase", "xray_debye", r_max=1e-5, q=p.q, intensity=p.intensity)


def test_powder_neutron_differs(dc):
    x = anatase_pattern(dc, radiation="xray")
    n = anatase_pattern(dc, radiation="neutron")
    np.testing.assert_allclose(x.q, n.q)
    assert not np.allclose(x.intensity / x.intensity.max(), n.intensity / n.intensity.max(),
                           atol=0.05)
    check_golden("api.powder_anatase", "neutron_debye", r_max=1e-5, q=n.q, intensity=n.intensity)


def test_powder_reloads_older_crystal(dc):
    """DISCUS holds one structure; powder() on an older Crystal reads it again."""
    big = dc.read_structure(ANATASE, supercell=(6, 6, 3))
    first = dc.powder(big, qmax=5.0)
    small = dc.read_structure(ANATASE, supercell=(2, 2, 1))
    assert dc.session().current_crystal is small
    again = dc.powder(big, qmax=5.0)
    assert dc.session().current_crystal is big
    np.testing.assert_allclose(again.intensity, first.intensity, rtol=1e-10)


def test_powder_bad_options(dc):
    cr = dc.read_structure(ANATASE)
    with pytest.raises(ValueError):
        dc.powder(cr, radiation="gamma")
    with pytest.raises(ValueError):
        dc.powder(cr, qmin=5, qmax=1)


# --------------------------------------------------------------- plotting
def test_plot(dc):
    mpl = pytest.importorskip("matplotlib")
    mpl.use("Agg")
    import matplotlib.pyplot as plt

    p = anatase_pattern(dc)
    ax = p.plot()
    line = ax.get_lines()[0]
    np.testing.assert_allclose(line.get_xdata(), p.q)
    ax2 = p.plot(x="tth")
    assert ax2.get_xlim()[0] == pytest.approx(p.two_theta()[0])
    plt.close("all")


@pytest.mark.parametrize("script", sorted(p.name for p in EXAMPLES.glob("*.py")))
def test_example_script(script, tmp_path):
    """Each script in python/examples/ runs on its own (own process, no display)."""
    pytest.importorskip("matplotlib")
    if not os.environ.get("DIFFUSE_PYTHON_PATH"):
        pytest.skip("DIFFUSE_PYTHON_PATH not set")
    env = os.environ | {"PYTHONPATH": str(REPO / "python"), "MPLBACKEND": "Agg"}
    r = subprocess.run([sys.executable, str(EXAMPLES / script)], cwd=tmp_path,
                       env=env, capture_output=True, text=True, timeout=300)
    assert r.returncode == 0, r.stdout + r.stderr
