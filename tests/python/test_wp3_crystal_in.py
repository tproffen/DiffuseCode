"""WP3, second slice: a crystal from Python arrays into DISCUS (no file),
and the ASE bridge.

The reference is always DISCUS itself: the same structure read from its
file. A crystal made from the arrays must hold the same data (round trip
Python -> Fortran -> Python) and give the same powder pattern, Fourier map
and PDF.

    DIFFUSE_PYTHON_PATH=<build>/python pytest tests/python/test_wp3_crystal_in.py
"""

from __future__ import annotations

import numpy as np
import pytest

from test_wp3_structures import ANATASE, STRUCTURES, T

NICKEL = T / "03.Pair_Distribution_Function/ni.cell"
CUO = T / "02.Calculation_Scattering/cuo.cell"


@pytest.fixture
def dc(suite):
    import diffuse
    return diffuse


def snapshot(dc):
    """Everything the C interface gives about the crystal DISCUS holds."""
    core = dc.session()._c
    pos, types, prop, mole = core.atoms()
    names, b_iso, occ = core.atom_types()
    return dict(cell=core.cell(), pos=pos, types=types, prop=prop, names=names,
                b_iso=b_iso, occ=occ)


def remake(dc, crystal):
    """make_crystal from the arrays of a crystal read from a file (B and
    occupancy per atom: two types may share a name)."""
    names, xyz = crystal.atoms()
    _pos, types, _prop, _mole = dc.session()._c.atoms()
    _n, b, occ = dc.session()._c.atom_types()
    return dc.make_crystal(crystal.lattice, names, xyz, space_group=crystal.space_group,
                           b_iso=b[types], occupancy=occ[types], cells=crystal.supercell,
                           types=types)


@pytest.mark.parametrize("path, spcgr, natoms", STRUCTURES, ids=lambda v: getattr(v, "stem", None))
def test_round_trip(dc, path, spcgr, natoms):
    """File -> arrays -> make_crystal -> arrays: identical data."""
    from_file = dc.read_structure(path, supercell=(2, 1, 2))
    ref = snapshot(dc)
    made = remake(dc, from_file)
    got = snapshot(dc)
    assert made.natoms == from_file.natoms and made.space_group == from_file.space_group
    assert got["cell"] == ref["cell"]
    np.testing.assert_array_equal(got["pos"], ref["pos"])
    np.testing.assert_array_equal(got["types"], ref["types"])
    np.testing.assert_array_equal(got["names"], ref["names"])
    np.testing.assert_array_equal(got["b_iso"], ref["b_iso"])
    np.testing.assert_array_equal(got["occ"], ref["occ"])


@pytest.mark.parametrize("path", [NICKEL, CUO], ids=lambda p: p.stem)
def test_same_calculations_as_from_file(dc, path):
    """Powder (Debye), Fourier plane and PDF of the crystal made from arrays
    equal those of the crystal read from the file."""
    from_file = dc.read_structure(path, supercell=(3, 3, 3))
    made = remake(dc, from_file)
    for crystal_a, crystal_b in ((from_file, made),):
        kw = dict(qmin=1.0, qmax=6.0, dq=0.01, radiation="xray")
        pa, pb = dc.powder(crystal_a, **kw), dc.powder(crystal_b, **kw)
        np.testing.assert_allclose(pb.intensity, pa.intensity, rtol=1e-6)
        kw = dict(ll=(0, 0, 0), lr=(4, 0, 0), ul=(0, 4, 0), npoints=(41, 41))
        fa, fb = dc.fourier(crystal_a, **kw), dc.fourier(crystal_b, **kw)
        np.testing.assert_allclose(fb.intensity, fa.intensity, rtol=1e-6, atol=1e-9 * fa.intensity.max())
        kw = dict(rmax=8.0, dr=0.02)
        ga, gb = dc.pdf(crystal_a, **kw), dc.pdf(crystal_b, **kw)
        np.testing.assert_allclose(gb.G, ga.G, rtol=1e-6, atol=1e-9)


def test_reloads_when_another_structure_was_read(dc):
    """A crystal made from arrays is sent again before a calculation if
    DISCUS holds another structure meanwhile."""
    xyz = np.array([[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]], float)
    cells = np.array([[i, j, k] for i in range(4) for j in range(4) for k in range(4)], float)
    xyz = (cells[:, None, :] + xyz[None]).reshape(-1, 3)
    cu = dc.make_crystal((3.61, 3.61, 3.61, 90, 90, 90), ["CU"] * len(xyz), xyz,
                         space_group="F m -3 m", b_iso=0.6, cells=(4, 4, 4))
    first = dc.powder(cu, qmin=2, qmax=5, dq=0.02)
    dc.read_structure(CUO, supercell=(2, 2, 2))
    again = dc.powder(cu, qmin=2, qmax=5, dq=0.02)
    np.testing.assert_allclose(again.intensity, first.intensity, rtol=1e-9)
    assert dc.session().current_crystal is cu
    assert cu.atom_types()[0].tolist() == ["CU"] and cu.atom_types()[1].tolist() == [0.6]


def test_types_in_order_of_first_appearance(dc):
    cr = dc.make_crystal((5, 5, 5, 90, 90, 90), ["O", "ti", "O", "Ti"],
                         [[0, 0, 0], [0.5, 0.5, 0.5], [0.5, 0, 0], [0, 0.5, 0.5]])
    names, xyz = cr.atoms()
    assert cr.atom_types()[0].tolist() == ["O", "TI"]
    assert names.tolist() == ["O", "TI", "O", "TI"]


def test_same_name_different_b_is_two_types(dc):
    """As in DISCUS: a type is name, B and occupancy."""
    cr = dc.make_crystal((5, 5, 5, 90, 90, 90), ["O", "O", "TI"],
                         [[0, 0, 0], [0.5, 0, 0], [0.5, 0.5, 0.5]], b_iso=[0.5, 0.9, 0.3])
    names, b, occ = cr.atom_types()
    assert names.tolist() == ["O", "O", "TI"] and b.tolist() == [0.5, 0.9, 0.3]


@pytest.mark.parametrize("kwargs, match", [
    (dict(names=["NI"] * 3), "shape"),
    (dict(names=["NICKEL"] * 4), "4 characters"),
    (dict(cells=(3, 1, 1)), "evenly"),
    (dict(b_iso={"CU": 1.0}), "no value"),
    (dict(b_iso=[1.0, 2.0]), "one per atom"),
])
def test_bad_arrays(dc, kwargs, match):
    args = dict(names=["NI"] * 4, cells=(1, 1, 1))
    args.update(kwargs)
    with pytest.raises(ValueError, match=match):
        dc.make_crystal((3.5, 3.5, 3.5, 90, 90, 90), args.pop("names"),
                        np.zeros((4, 3)), **args)


def test_unknown_space_group_is_a_discus_error(dc):
    with pytest.raises(dc.MacroError, match="Space group"):
        dc.make_crystal((3.5, 3.5, 3.5, 90, 90, 90), ["NI"], [[0, 0, 0]], space_group="Q 99")
    # DISCUS is still usable afterwards
    assert dc.read_structure(NICKEL).natoms == 4


# ------------------------------------------------------------------ ASE
ase = pytest.importorskip("ase")


def test_from_ase_bulk(dc):
    from ase.build import bulk
    atoms = bulk("Cu", "fcc", a=3.61, cubic=True) * (2, 3, 1)
    cr = dc.from_ase(atoms, b_iso=0.6)
    assert cr.natoms == 24 and cr.space_group == "P 1"     # DISCUS's spelling
    assert cr.lattice == pytest.approx((7.22, 10.83, 3.61, 90, 90, 90))
    np.testing.assert_allclose(cr.cartesian(), atoms.get_positions(), atol=1e-9)
    assert set(cr.atoms()[0]) == {"CU"}


def test_to_ase_and_back(dc):
    """read -> to_ase -> from_ase: same atoms (the ASE cell is the whole
    supercell, so the new crystal is one cell of that size) and the same
    Debye powder pattern, which depends only on the distances."""
    cr = dc.read_structure(ANATASE, supercell=(3, 3, 2))
    atoms = cr.to_ase()
    assert len(atoms) == cr.natoms
    assert sorted(set(atoms.get_chemical_symbols())) == ["O", "Ti"]
    assert list(atoms.arrays["discus_names"][:1]) == [cr.atoms()[0][0]]
    np.testing.assert_allclose(atoms.cell.cellpar()[:3], np.array(cr.lattice[:3]) * (3, 3, 2))
    p1 = dc.powder(cr, qmin=1, qmax=6, dq=0.02)
    back = dc.from_ase(atoms, b_iso=dict(zip(*cr.atom_types()[:2])))
    assert back.natoms == cr.natoms and back.atoms()[0].tolist() == cr.atoms()[0].tolist()
    np.testing.assert_allclose(back.cartesian(), cr.cartesian(), atol=1e-8)
    p2 = dc.powder(back, qmin=1, qmax=6, dq=0.02)
    np.testing.assert_allclose(p2.intensity, p1.intensity, rtol=1e-6)


@pytest.mark.filterwarnings("ignore:crystal system:UserWarning")   # ASE's note on the CIF setting
def test_ase_reads_the_same_cif_as_discus(dc):
    """ASE's CIF reader and DISCUS's give the same anatase unit cell: an
    independent check of the space-group expansion, through from_ase."""
    import ase.io
    atoms = ase.io.read(ANATASE)
    via_ase = dc.from_ase(atoms)
    sites_ase = sorted((n, *np.round(np.mod(x, 1), 4)) for n, x in zip(*via_ase.atoms()))
    discus = dc.read_structure(ANATASE)
    sites_discus = sorted((n, *np.round(np.mod(x, 1), 4)) for n, x in zip(*discus.atoms()))
    elements = {"TI4+": "TI", "O2-": "O"}
    assert [(elements.get(s[0], s[0]),) + s[1:] for s in sites_discus] == sites_ase
