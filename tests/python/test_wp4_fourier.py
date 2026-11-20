"""WP4 safety net: single-crystal Fourier options, checked against physics.

Golden data for the WP4 work package (the Fourier shim). Where possible the
reference is physics rather than an earlier DISCUS run:

* fcc extinction rules and equal intensities (nickel, neutrons, no ADPs);
* the Debye-Waller factor exp(-B |d*|^2 / 2);
* X-ray form factors against the International Tables (IT92, via gemmi);
* anomalous dispersion: Friedel pairs differ in non-centrosymmetric ZnSe by
  the amount the structure factor predicts, and stay equal in nickel;
* energy and wavelength settings agree; seeded 'lots' runs repeat exactly;
* the NUFFT technique against the standard one.

The DISCUS commands are run directly (diffuse.run_macro), so options that
the Python API does not offer yet are covered too.

    DIFFUSE_PYTHON_PATH=<build>/python pytest tests/python/test_wp4_fourier.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from compare import check_golden, rvalue

REPO = Path(__file__).resolve().parents[2]
NI = REPO / "tests/03.Pair_Distribution_Function/ni.cell"                 # Fm-3m, B = 0.2
ZNSE = REPO / "tests/08.Nanoparticle_example/COMPLEX/CELL/znse_zincblende.cell"  # F-43m
LAYER = REPO / "tests/02.Calculation_Scattering/example3.stru"           # 50x50 Cu/Au layer
A_NI, B_NI = 3.5234, 0.2
HK0 = dict(ll="-3,-3,0", lr="3,-3,0", ul="-3,3,0", na=61, no=61)


@pytest.fixture
def dc(suite):
    import diffuse
    return diffuse


def fourier(dc, *, ll, lr, ul=None, na, no=1, radiation="neutron", temp="ign",
            extra=""):
    """Run the fourier menu with the given settings; return the KUPLOT data.

    A line gives (x, intensity); a plane gives (x, y, intensity[na, no]).
    ``extra`` holds further fourier-menu commands (wvle, disp, delf, set ...).
    """
    s = dc.session()
    out = s.workdir / "wp4.inte"
    out.unlink(missing_ok=True)
    s.run(f"""
discus
  fourier
    reset
    {radiation}
    ll {ll}
    lr {lr}
    ul {ul or ll}
    na {na}
    no {no}
    abs h
    ord k
    temp {temp}
    disp off
    set aver, 0
    lots off
{extra}
    run
  exit
  output
    outf {out.name}
    value inte
    form stan
    run
  exit
exit
kuplot
  rese
  load {"xy" if no == 1 else "ni"}, {out.name}
exit
""", cwd=s.workdir)
    if no == 1:
        x, i = s.dataset(1)
        return x.astype(float), i.astype(float)
    x, y, z = s.dataset_2d(1)
    return x.astype(float), y.astype(float), z.astype(float)


# ------------------------------------------------------------ basic physics
def test_fcc_extinction_and_equal_intensities(dc):
    """Ni, 4x4x4 cells, neutrons, no ADPs: hk0 reflections only for h, k even,
    all of the same intensity; zero at mixed parity and between Bragg peaks."""
    dc.read_structure(NI, supercell=4)
    h, k, z = fourier(dc, ll="0,0,0", lr="4,0,0", ul="0,4,0", na=9, no=9)
    hh, kk = np.meshgrid(h, k, indexing="ij")
    allowed = (np.isclose(hh % 2, 0) & np.isclose(kk % 2, 0))
    np.testing.assert_allclose(z[allowed], z.max(), rtol=1e-5)
    assert np.all(z[~allowed] < 1e-6 * z.max())


def test_debye_waller_factor(dc):
    """With ADPs, each reflection is scaled by exp(-B |d*|^2 / 2)."""
    dc.read_structure(NI, supercell=4)
    plane = dict(ll="0,0,0", lr="4,0,0", ul="0,4,0", na=9, no=9)
    h, k, cold = fourier(dc, temp="ign", **plane)
    _, _, warm = fourier(dc, temp="use", **plane)
    hh, kk = np.meshgrid(h, k, indexing="ij")
    peak = cold > 1e-6 * cold.max()
    d2 = (hh**2 + kk**2) / A_NI**2
    np.testing.assert_allclose(warm[peak] / cold[peak], np.exp(-B_NI * d2[peak] / 2), rtol=1e-3)


def test_xray_form_factor_matches_it92(dc):
    """X-ray / neutron intensity follows f(s)^2 of the International Tables."""
    gemmi = pytest.importorskip("gemmi")
    dc.read_structure(NI, supercell=4)
    plane = dict(ll="0,0,0", lr="4,0,0", ul="0,4,0", na=9, no=9)
    h, k, neutron = fourier(dc, radiation="neutron", **plane)
    _, _, xray = fourier(dc, radiation="xray", **plane)
    hh, kk = np.meshgrid(h, k, indexing="ij")
    peak = neutron > 1e-6 * neutron.max()
    s = np.sqrt((hh**2 + kk**2)[peak]) / A_NI / 2              # sin(theta)/lambda
    f = np.array([gemmi.Element("Ni").it92.calculate_sf(float(v) ** 2) for v in s])
    ratio = xray[peak] / neutron[peak] / f**2
    assert np.ptp(ratio) / ratio.mean() < 0.01


def test_energy_equals_wavelength(dc):
    """'energy E' (keV) gives the same map as 'wvle 12.3984/E'."""
    dc.read_structure(ZNSE, supercell=3)
    *_, a = fourier(dc, radiation="xray", extra="    wvle 1.0", **HK0)
    *_, b = fourier(dc, radiation="xray", extra="    energy 12.398419843", **HK0)
    assert rvalue(a, b) < 1e-6


# -------------------------------------------------------- anomalous dispersion
ANOMALOUS = "    wvle 1.0\n    delf ZN, -1.5, 0.7\n    delf SE, -3.0, 0.5\n    disp anom"
LINE_111 = dict(ll="-3,-3,-3", lr="3,3,3", na=61)          # h = k = l from -3 to 3


def test_without_dispersion_friedel_pairs_are_equal(dc):
    dc.read_structure(ZNSE, supercell=3)
    _, i = fourier(dc, radiation="xray", extra="    wvle 1.0", **LINE_111)
    np.testing.assert_allclose(i, i[::-1], rtol=1e-6)


def test_dispersion_breaks_friedel_law_in_znse(dc):
    """Zincblende: F(hhh) = 4 (f_Zn + f_Se exp(i pi 3h/2)). With f = f0 + f' + i f''
    the 111 and -1-1-1 intensities differ by a predictable ratio."""
    gemmi = pytest.importorskip("gemmi")
    dc.read_structure(ZNSE, supercell=3)
    x, i = fourier(dc, radiation="xray", extra=ANOMALOUS, **LINE_111)
    i111, i_111 = i[np.argmin(abs(x - 1))], i[np.argmin(abs(x + 1))]
    s = np.sqrt(3) / 5.496 / 2
    f_zn = gemmi.Element("Zn").it92.calculate_sf(s**2) - 1.5 + 0.7j
    f_se = gemmi.Element("Se").it92.calculate_sf(s**2) - 3.0 + 0.5j
    expected = abs(f_zn - 1j * f_se) ** 2 / abs(f_zn + 1j * f_se) ** 2
    assert i111 / i_111 == pytest.approx(expected, rel=5e-3)
    assert abs(i111 / i_111 - 1) > 0.01


def test_dispersion_keeps_friedel_law_in_centrosymmetric_nickel(dc):
    dc.read_structure(NI, supercell=3)
    _, i = fourier(dc, radiation="xray",
                   extra="    wvle 1.0\n    delf NI, -2.0, 1.0\n    disp anom", **LINE_111)
    assert rvalue(i, i[::-1]) < 1e-4


def test_disp_on_is_not_anomalous(dc):
    """Only 'disp anom' switches dispersion on; 'disp on' leaves it off."""
    dc.read_structure(ZNSE, supercell=3)
    _, i = fourier(dc, radiation="xray",
                   extra=ANOMALOUS.replace("disp anom", "disp on"), **LINE_111)
    np.testing.assert_allclose(i, i[::-1], rtol=1e-6)


def test_internal_dispersion_tables_golden(dc):
    """Tabulated f', f'' for Mo K-alpha1 ('wvle moa1', 'delf all, internal')."""
    dc.read_structure(ZNSE, supercell=3)
    x, i = fourier(dc, radiation="xray",
                   extra="    wvle moa1\n    delf all, internal\n    disp anom", **LINE_111)
    assert rvalue(i, i[::-1]) > 0.005                     # Friedel pairs differ
    check_golden("wp4.fourier", "znse_111_moa1_anomalous", r_max=1e-5, h=x, intensity=i)


# ------------------------------------------------------- radiation goldens
@pytest.mark.parametrize("radiation", ["xray", "neutron", "electron"])
def test_radiation_golden(dc, radiation):
    dc.read_structure(ZNSE, supercell=3)
    h, k, z = fourier(dc, radiation=radiation, temp="use", extra="    wvle 0.71", **HK0)
    assert np.all(z >= 0)
    check_golden("wp4.fourier", f"znse_hk0_{radiation}", r_max=1e-5, h=h, k=k, intensity=z)


# --------------------------------------------------------------- technique
def test_nufft_matches_standard_for_layer(dc):
    """NUFFT and the standard technique agree for the 50x50 layer on the grid
    that 04.SRO/SRO_3 uses."""
    dc.read_structure(LAYER)
    plane = dict(ll="-10,-10,0", lr="10,-10,0", ul="-10,10,0", na=201, no=201)
    *_, standard = fourier(dc, radiation="xray", extra="    set technique:turbo", **plane)
    *_, nufft = fourier(dc, radiation="xray", extra="    set technique:nufft", **plane)
    assert rvalue(standard, nufft) < 1e-4


@pytest.mark.xfail(strict=True, reason=(
    "NUFFT gives wrong maps for small crystals on many grids: Ni 4x4x4, hk0 "
    "-3..3 with 61x61 points has R = 1.7 against the standard technique, with "
    "spurious maxima every 0.08-0.16 in h. Some size/grid combinations are "
    "exact (10x10 cells, or a 9x9 grid over 0..4). Remove the mark once fixed."))
def test_nufft_matches_standard_for_small_crystal(dc):
    dc.read_structure(NI, supercell=4)
    *_, standard = fourier(dc, extra="    set technique:turbo", **HK0)
    *_, nufft = fourier(dc, extra="    set technique:nufft", **HK0)
    assert rvalue(standard, nufft) < 1e-4


# -------------------------------------------------------------------- lots
LOTS = "    lots box, 10, 10, 1, 20, yes"
LOTS_PLANE = dict(ll="0,0,0", lr="2,0,0", ul="0,2,0", na=81, no=81)


def lots_map(dc, seed):
    dc.read_structure(LAYER)
    dc.run_macro(f"seed {seed}")
    return fourier(dc, radiation="xray", extra=LOTS, **LOTS_PLANE)


def test_lots_with_seed_repeat_exactly(dc):
    *_, first = lots_map(dc, 4321)
    *_, again = lots_map(dc, 4321)
    *_, other = lots_map(dc, 999)
    assert np.all(np.isfinite(first))
    np.testing.assert_array_equal(first, again)
    assert rvalue(first, other) > 0.01                     # another seed, other lots


def test_lots_golden(dc):
    h, k, z = lots_map(dc, 4321)
    # Differs between compilers/CPUs by R ~ 4e-4 (gfortran 11 vs 13 with
    # -ffast-math -march=native); same build: identical (test above).
    check_golden("wp4.fourier", "layer_lots_seed4321", r_max=2e-3, h=h, k=k, intensity=z)
