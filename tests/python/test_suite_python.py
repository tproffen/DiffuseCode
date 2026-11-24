"""Tests for the in-process Python interface (the session that runs DISCUS;
named after the f2py module suite_python that came first).

These run DISCUS inside the Python process and read results back as NumPy
arrays. Each calculation is checked two ways:

* against the legacy reference file in ``tests/<folder>/RESULTS/`` (the same
  file the macro test compares with), which shows that the Python path gives
  the same answer as ``discus_suite``;
* against golden data in ``tests/golden/<id>/`` (see compare.py), which is
  what the Python API work (WP3–WP5) will be tested against.

    DIFFUSE_PYTHON_PATH=<build>/python pytest tests/python/test_suite_python.py

Without the module these tests are skipped (set DIFFUSE_REQUIRE_PYTHON=1 to
make that an error instead, as CI does).
"""

from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pytest

from compare import check_golden, load_columns, rvalue
from diffuse import MacroError

TESTS = Path(__file__).resolve().parents[1]
SCATTERING = TESTS / "02.Calculation_Scattering"


@pytest.fixture
def scattering_dir(tmp_path):
    """Fresh folder with the inputs from 02.Calculation_Scattering."""
    for name in ("silicon.cell", "primitive.cell", "example3.stru"):
        shutil.copy(SCATTERING / name, tmp_path)
    shutil.copytree(SCATTERING / "RESULTS", tmp_path / "RESULTS")
    return tmp_path


# ------------------------------------------------------------------ basics
def test_expression(suite, tmp_path):
    assert suite.value("2*3+1") == 7.0
    suite.run("r[1] = 2.5\nr[2] = r[1]**2", tmp_path)
    assert suite.value("r[2]") == pytest.approx(6.25)


def test_macro_error_raises(suite, tmp_path):
    with pytest.raises(MacroError):
        suite.run("discus\n  no_such_command 1,2\nexit", tmp_path)
    # the session is still usable afterwards
    suite.run("r[3] = 4", tmp_path)
    assert suite.value("r[3]") == 4.0


def test_missing_dataset_raises(suite, tmp_path):
    suite.run("kuplot\n  rese\nexit", tmp_path)
    with pytest.raises(IndexError):
        suite.dataset(1)


# --------------------------------------------------------------- structure
def test_build_crystal(suite, scattering_dir):
    suite.run("""
discus
  read
    cell silicon.cell, 3, 2, 1
  r[10] = n[1]
  r[11] = lat[1]
  r[12] = lat[6]
exit
""", scattering_dir)
    assert suite.value("r[10]") == 6            # one atom per cell, 3*2*1 cells
    assert suite.value("r[11]") == pytest.approx(5.0)
    assert suite.value("r[12]") == pytest.approx(90.0)


# ----------------------------------------------------------------- Fourier
FOURIER_1D = """
discus
  read
    cell silicon.cell, 7, 1, 1
  fourier
    reset
    neutron
    ll -2.5, 0.00, 0.00
    lr  2.5, 0.00, 0.00
    ul -2.5, 0.00, 0.00
    na 1001
    no 1
    abs h
    ord k
    temp ign
    disp off
    set aver, 0.0
    lots off
    run
  exit
  output
    outf py_1d.inte
    value inte
    form stan
    run
  exit
exit
kuplot
  rese
  load xy, py_1d.inte
exit
"""


def test_fourier_1d(suite, scattering_dir):
    """Same calculation as tests/02.Calculation_Scattering/1d.mac 7."""
    suite.run(FOURIER_1D, scattering_dir)
    h, inte = suite.dataset(1)
    assert h.shape == (1001,)
    assert h[0] == pytest.approx(-2.5) and h[-1] == pytest.approx(2.5)
    # a row of atoms along a: the strongest intensity is at an integer h
    assert h[np.argmax(inte)] == pytest.approx(round(h[np.argmax(inte)]), abs=1e-4)

    legacy = load_columns(scattering_dir / "RESULTS" / "1d.inte")
    np.testing.assert_allclose(h, legacy[:, 0], atol=1e-5)
    assert rvalue(legacy[:, 1], inte) < 1e-3

    check_golden("02.1d", "fourier_1d", r_max=1e-5, h=h, intensity=inte)


def test_fourier_2d(suite, scattering_dir):
    """2-D (hk0) section of example3.stru, coarser than 2d.mac (121 x 121)."""
    suite.run("""
discus
  read
    stru example3.stru
  fourier
    reset
    xray
    ll 0.0, 0.0, 0.0
    lr 4.0, 0.0, 0.0
    ul 0.0, 4.0, 0.0
    na 121
    no 121
    temp use
    disp off
    wvle 1.00
    abs h
    ord k
    set aver, 0
    run
  exit
  output
    outf py_2d.inte
    value inte
    form stan
    run
  exit
exit
kuplot
  rese
  load ni, py_2d.inte
exit
""", scattering_dir)
    h, k, inte = suite.dataset_2d(1)
    assert h.shape == (121,) and k.shape == (121,) and inte.shape == (121, 121)
    assert h[-1] == pytest.approx(4.0) and k[-1] == pytest.approx(4.0)
    # strongest intensity at the origin and at the Bragg points (integer h,k)
    i, j = np.unravel_index(np.argmax(inte), inte.shape)
    assert h[i] == pytest.approx(round(h[i]), abs=1e-4)
    assert k[j] == pytest.approx(round(k[j]), abs=1e-4)

    check_golden("02.2d", "fourier_2d", r_max=1e-5, h=h, k=k, intensity=inte)


# ------------------------------------------------------------------ powder
def test_powder_debye(suite, scattering_dir):
    """Same calculation as tests/02.Calculation_Scattering/powder.mac 10."""
    suite.run("""
set parallel, useomp:serial
discus
  read
    cell primitive.cell, 10, 10, 10
  powder
    reset
    xray
    set axis, q
    set calc, debye
    set disp, off
    set delta, 0.0
    set qmin, 0.500
    set qmax, 7.100
    set dq, 0.005
    set profile, off
    set temp, use
    set wvle, 0.20000
    set four, four
    set lpcor, bragg, 26.58
    run
  exit
  output
    outf py_powder.inte
    value inte
    form powder, tth, 1.0, 180, 0.01
    run
  exit
exit
kuplot
  rese
  load xy, py_powder.inte
exit
""", scattering_dir)
    tth, inte = suite.dataset(1)
    assert np.all(np.diff(tth) > 0)
    assert np.all(inte >= 0)

    legacy = load_columns(scattering_dir / "RESULTS" / "powder.inte")
    np.testing.assert_allclose(tth, legacy[:, 0], atol=1e-4)
    assert rvalue(legacy[:, 1], inte) < 1e-3

    check_golden("02.powder", "powder_debye", r_max=1e-5, tth=tth, intensity=inte)
