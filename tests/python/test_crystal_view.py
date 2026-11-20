"""Crystal geometry (lattice vectors, Cartesian positions) and Crystal.view().

view() is checked with a stand-in for py3Dmol that records the calls (what
is drawn), plus one smoke test with the real py3Dmol when it is installed.

    DIFFUSE_PYTHON_PATH=<build>/python pytest tests/python/test_crystal_view.py
"""

from __future__ import annotations

import sys
import types
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[2]
EXAMPLES = REPO / "python" / "examples"
HEX = REPO / "tests" / "03.Pair_Distribution_Function" / "water_hex.cell"


@pytest.fixture
def dc(suite):
    import diffuse
    return diffuse


# ---------------------------------------------------------------- geometry
def test_lattice_vectors_hexagonal(dc):
    cr = dc.read_structure(HEX)
    a, b, c, al, be, ga = cr.lattice
    L = cr.lattice_vectors()
    np.testing.assert_allclose(np.linalg.norm(L, axis=1), (a, b, c), rtol=1e-12)
    cos = lambda u, v: u @ v / np.linalg.norm(u) / np.linalg.norm(v)
    np.testing.assert_allclose(np.degrees(np.arccos([cos(L[1], L[2]), cos(L[0], L[2]),
                                                     cos(L[0], L[1])])), (al, be, ga))
    assert np.linalg.det(L) == pytest.approx(cr.volume, rel=1e-10)
    assert L[0, 1] == L[0, 2] == L[1, 2] == 0           # a along x, b in xy


def test_cartesian_nearest_neighbours_nickel(dc):
    cr = dc.read_structure(EXAMPLES / "nickel.cell", supercell=3)
    xyz = cr.cartesian()
    assert xyz.shape == (cr.natoms, 3)
    d = np.linalg.norm(xyz[1:] - xyz[0], axis=1)
    assert d.min() == pytest.approx(cr.a / np.sqrt(2), rel=1e-6)    # fcc
    assert np.count_nonzero(np.isclose(d, d.min())) <= 12


def test_element_names():
    from diffuse.crystal import _element
    assert [_element(n) for n in ["TI4+", "O2-", "CU", "NI", "D", "1"]] == \
        ["Ti", "O", "Cu", "Ni", "D", "X"]


# -------------------------------------------------------------------- view
class FakeView:
    """Records what Crystal.view() asks py3Dmol to draw."""

    def __init__(self, width, height):
        self.size = (width, height)
        self.models, self.styles, self.lines = [], [], []

    def addModel(self, data, fmt):
        self.models.append((data, fmt))

    def setStyle(self, sel, style):
        self.styles.append((sel, style))

    def addLine(self, spec):
        self.lines.append(spec)

    def setBackgroundColor(self, color):
        self.background = color

    def zoomTo(self):
        pass


@pytest.fixture
def fake3dmol(monkeypatch):
    monkeypatch.setitem(sys.modules, "py3Dmol", types.SimpleNamespace(view=FakeView))


def model_atoms(view):
    data, fmt = view.models[0]
    assert fmt == "xyz"
    lines = data.splitlines()
    assert int(lines[0]) == len(lines) - 2
    return [line.split() for line in lines[2:]]


def test_view_model_and_crystal_box(dc, fake3dmol):
    cr = dc.read_structure(EXAMPLES / "nickel.cell", supercell=3)
    v = cr.view(width=300, height=200)
    atoms = model_atoms(v)
    assert len(atoms) == cr.natoms and {a[0] for a in atoms} == {"Ni"}
    np.testing.assert_allclose(np.array([a[1:] for a in atoms], float), cr.cartesian(), atol=1e-4)
    assert v.size == (300, 200)
    # 12 edges of the 3 x 3 x 3 box, each 3 a long
    assert len(v.lines) == 12
    for line in v.lines:
        p, q = (np.array([line[k][c] for c in "xyz"]) for k in ("start", "end"))
        assert np.linalg.norm(q - p) == pytest.approx(3 * cr.a, rel=1e-6)


def test_view_box_of_a_stru_file(dc, fake3dmol):
    """50 x 50 x 1 cells of 5 A: the box spans the whole layer."""
    cr = dc.read_structure(EXAMPLES / "cuau_layer.stru")
    v = cr.view()
    lengths = sorted({round(float(np.linalg.norm([l["end"][c] - l["start"][c] for c in "xyz"])), 3)
                      for l in v.lines})
    assert lengths == [5.0, 250.0]


def test_view_options(dc, fake3dmol):
    cr = dc.read_structure(EXAMPLES / "cuau_layer.stru")
    names, _ = cr.atoms()
    v = cr.view(show="au", colors={"AU": "gold"}, cell="unit", scale=2, style="ballstick")
    atoms = model_atoms(v)
    assert len(atoms) == np.count_nonzero(names == "AU") and {a[0] for a in atoms} == {"Au"}
    base, colored = v.styles
    assert base == ({}, {"sphere": {"scale": 0.5}, "stick": {"radius": 0.24}})
    assert colored[1]["sphere"]["color"] == "gold"
    assert len(colored[0]["index"]) == len(atoms)
    assert len(v.lines) == 12 and all(            # one 5 A unit cell
        np.linalg.norm([l["end"][c] - l["start"][c] for c in "xyz"]) == pytest.approx(5.0)
        for l in v.lines)
    assert cr.view(cell=None).lines == []


def test_view_bad_options(dc, fake3dmol):
    cr = dc.read_structure(EXAMPLES / "nickel.cell")
    for kwargs in (dict(style="cartoon"), dict(cell="box"), dict(show="FE"), dict(scale=0)):
        with pytest.raises(ValueError):
            cr.view(**kwargs)


def test_view_with_py3dmol(dc):
    py3Dmol = pytest.importorskip("py3Dmol")
    cr = dc.read_structure(EXAMPLES / "anatase.cif", supercell=2)
    v = cr.view(style="stick")
    assert isinstance(v.viewer, py3Dmol.view)
    html = v.write_html()
    assert "addModel" in html and "addLine" in html
    assert v.addLabel("x", {"position": {"x": 0, "y": 0, "z": 0}}) is v.viewer   # passed on
