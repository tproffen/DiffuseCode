"""Tests for diffuse.fourier and diffuse.pdf.

Both are checked against independent references where they exist:

* PDF: the legacy reference files in tests/03.Pair_Distribution_Function
  (RESULTS/*.grcalc, made by pdf.mac from water_hex.cell);
* Fourier: tests/02.Calculation_Scattering/RESULTS/1d.inte and the golden
  2-D map of test_suite_python.py (02.2d).

    DIFFUSE_PYTHON_PATH=<build>/python pytest tests/python/test_fourier_pdf.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from compare import check_golden, load_columns, rvalue

REPO = Path(__file__).resolve().parents[2]
PDF_DIR = REPO / "tests" / "03.Pair_Distribution_Function"
SCAT_DIR = REPO / "tests" / "02.Calculation_Scattering"
GOLDEN = REPO / "tests" / "golden"


@pytest.fixture
def dc(suite):
    import diffuse
    return diffuse


# ------------------------------------------------------------------- PDF
@pytest.mark.parametrize("kwargs, rmax, dr, reference", [
    (dict(thermal=False), 10.5, 0.001, "no_corr.grcalc"),
    (dict(), 10.5, 0.001, "thermal.grcalc"),
    (dict(qmax=30.0), 10.5, 0.001, "thermal_qmax30.grcalc"),
    (dict(qmax=30.0, qdamp=0.02), 110.0, 0.01, "thermal_damp_020.grcalc"),
])
def test_pdf_matches_legacy_reference(dc, kwargs, rmax, dr, reference):
    """Same settings as tests/03.Pair_Distribution_Function/pdf.mac."""
    cr = dc.read_structure(PDF_DIR / "water_hex.cell")
    g = dc.pdf(cr, rmax=rmax, dr=dr, **kwargs)
    ref = load_columns(PDF_DIR / "RESULTS" / reference)
    np.testing.assert_allclose(g.r, ref[:, 0], atol=1e-6)
    assert rvalue(ref[:, 1], g.G) < 1e-6


def test_pdf_nickel_peaks_and_golden(dc):
    cr = dc.read_structure(PDF_DIR / "ni.cell")
    g = dc.pdf(cr, rmax=10.0, dr=0.01)
    assert g.r[0] == pytest.approx(0.01) and g.r[-1] == pytest.approx(10.0)
    # fcc nickel: first neighbours at a/sqrt(2)
    first = g.r[np.argmax(g.G * (g.r < 3.0))]
    assert first == pytest.approx(3.5234 / np.sqrt(2), abs=0.02)
    # below the first neighbour distance G(r) is the -4 pi rho0 r baseline
    rho0 = 4 / 3.5234**3
    low = g.r < 2.0
    np.testing.assert_allclose(g.G[low], -4 * np.pi * rho0 * g.r[low], atol=1e-3)
    check_golden("api.pdf_ni", "thermal", r_max=1e-6, r=g.r, G=g.G)


def test_pdf_neutron_and_xray_differ(dc):
    """Ti and O weigh differently for X-rays and neutrons (a one-element
    structure such as water_hex.cell, O only, would give the same G(r))."""
    cr = dc.read_structure(REPO / "python" / "examples" / "anatase.cif")
    x = dc.pdf(cr, rmax=8, dr=0.01, radiation="xray")
    n = dc.pdf(cr, rmax=8, dr=0.01, radiation="neutron")
    assert rvalue(x.G, n.G) > 0.05


def test_pdf_finite_particle(dc):
    """With particle_diameter, G(r) vanishes beyond the particle size."""
    a = 3.5234
    cr = dc.read_structure(PDF_DIR / "ni.cell", supercell=5)      # 5a cube
    nano = dc.pdf(cr, rmax=40, dr=0.02, particle_diameter=5 * a)
    bulk = dc.pdf(cr, rmax=40, dr=0.02)
    far = nano.r > 5 * a * np.sqrt(3) + 1        # beyond the cube diagonal
    assert np.abs(nano.G[far]).max() < 1e-9
    assert np.abs(bulk.G[far]).max() > 1.0


def test_pdf_bad_options(dc):
    cr = dc.read_structure(PDF_DIR / "ni.cell")
    with pytest.raises(ValueError):
        dc.pdf(cr, radiation="muon")
    with pytest.raises(ValueError):
        dc.pdf(cr, rmax=1.0, dr=2.0)
    with pytest.raises(ValueError):
        dc.pdf(cr, particle_diameter=-1)


# ---------------------------------------------------------------- Fourier
def test_fourier_line_matches_legacy(dc):
    """Same calculation as tests/02.Calculation_Scattering/1d.mac 7."""
    cr = dc.read_structure(SCAT_DIR / "silicon.cell", supercell=(7, 1, 1))
    f = dc.fourier(cr, ll=(-2.5, 0, 0), lr=(2.5, 0, 0), npoints=1001,
                   radiation="neutron", temperature=False)
    assert f.is_line and f.intensity.shape == (1001,) and f.x_label == "h"
    ref = load_columns(SCAT_DIR / "RESULTS" / "1d.inte")
    np.testing.assert_allclose(f.x, ref[:, 0], atol=1e-5)
    assert rvalue(ref[:, 1], f.intensity) < 1e-5
    np.testing.assert_allclose(f.hkl()[[0, -1]], [[-2.5, 0, 0], [2.5, 0, 0]])


def test_fourier_plane_matches_golden(dc):
    """Same map as test_suite_python.test_fourier_2d (golden 02.2d)."""
    cr = dc.read_structure(SCAT_DIR / "example3.stru")
    f = dc.fourier(cr, ll=(0, 0, 0), lr=(4, 0, 0), ul=(0, 4, 0), npoints=(121, 121),
                   wavelength=1.0)
    ref = np.load(GOLDEN / "02.2d" / "fourier_2d.npz")
    np.testing.assert_allclose(f.x, ref["h"], atol=1e-6)
    np.testing.assert_allclose(f.y, ref["k"], atol=1e-6)
    assert rvalue(ref["intensity"], f.intensity) < 1e-5


def test_fourier_orientation(dc):
    """intensity[i, j]: i runs along ll->lr, j along ll->ul, and hkl() agrees.

    A row of atoms along a, neutrons, no displacement parameters: the
    intensity varies with h only, and is exactly constant along k.
    """
    cr = dc.read_structure(SCAT_DIR / "silicon.cell", supercell=(7, 1, 1))
    f = dc.fourier(cr, ll=(-1.5, -1, 0), lr=(1.5, -1, 0), ul=(-1.5, 1, 0),
                   npoints=(61, 41), radiation="neutron", temperature=False)
    assert f.intensity.shape == (61, 41)
    assert (f.x_label, f.y_label) == ("h", "k")
    assert np.ptp(f.intensity[:, 0]) > 0.5 * f.intensity.max()
    assert np.ptp(f.intensity[30, :]) < 1e-6 * f.intensity.max()
    hkl = f.hkl()
    assert hkl.shape == (61, 41, 3)
    np.testing.assert_allclose(hkl[0, 0], (-1.5, -1, 0))
    np.testing.assert_allclose(hkl[-1, 0], (1.5, -1, 0))
    np.testing.assert_allclose(hkl[0, -1], (-1.5, 1, 0))
    i, _ = np.unravel_index(np.argmax(f.intensity), f.intensity.shape)
    assert hkl[i, 0, 0] == pytest.approx(round(hkl[i, 0, 0]), abs=1e-6)


def test_fourier_diagonal_plane_axes(dc):
    cr = dc.read_structure(SCAT_DIR / "example3.stru")
    f = dc.fourier(cr, ll=(0, 0, 0), lr=(3, 3, 0), ul=(0, 0, 1), npoints=(31, 11))
    assert (f.x_label, f.y_label) == ("h", "l")
    np.testing.assert_allclose(f.hkl()[-1, -1], (3, 3, 1))


def test_fourier_subtract_average_leaves_diffuse(dc):
    """Removing the average structure takes out the Bragg peaks of the Cu/Au layer."""
    cr = dc.read_structure(SCAT_DIR / "example3.stru")
    plane = dict(ll=(0, 0, 0), lr=(2, 0, 0), ul=(0, 2, 0), npoints=(81, 81))
    total = dc.fourier(cr, **plane)
    diffuse = dc.fourier(cr, subtract_average=100, **plane)
    assert diffuse.intensity.max() < 0.01 * total.intensity.max()
    assert diffuse.intensity.mean() > 0
    check_golden("api.fourier_cuau", "diffuse_hk0", r_max=1e-5,
                 h=diffuse.x, k=diffuse.y, intensity=diffuse.intensity)


def test_fourier_bad_options(dc):
    cr = dc.read_structure(SCAT_DIR / "silicon.cell")
    with pytest.raises(ValueError):
        dc.fourier(cr, ll=(0, 0, 0), lr=(0, 0, 0))
    with pytest.raises(ValueError):
        dc.fourier(cr, ll=(0, 0), lr=(1, 0))
    with pytest.raises(ValueError):
        dc.fourier(cr, ll=(0, 0, 0), lr=(1, 0, 0), radiation="gamma")
    with pytest.raises(ValueError):
        dc.fourier(cr, ll=(0, 0, 0), lr=(1, 0, 0), subtract_average=150)


# --------------------------------------------------------------- plotting
def test_plots(dc):
    mpl = pytest.importorskip("matplotlib")
    mpl.use("Agg")
    import matplotlib.pyplot as plt

    cr = dc.read_structure(SCAT_DIR / "silicon.cell", supercell=(7, 1, 1))
    line = dc.fourier(cr, ll=(-2, 0, 0), lr=(2, 0, 0), npoints=101)
    assert len(line.plot(log=True).get_lines()) == 1
    plane = dc.fourier(cr, ll=(-2, -2, 0), lr=(2, -2, 0), ul=(-2, 2, 0), npoints=21)
    ax = plane.plot(log=True)
    assert len(ax.images) == 1 and ax.images[0].get_array().shape == (21, 21)
    g = dc.pdf(dc.read_structure(PDF_DIR / "ni.cell"), rmax=5)
    assert len(g.plot().get_lines()) >= 1
    plt.close("all")
