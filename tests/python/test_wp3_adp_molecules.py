"""WP3, third slice: anisotropic displacement parameters and molecules,
out of DISCUS and into it.

References: DISCUS reading the same data from a file (round trips must be
identical, and Fourier maps, which depend on the ADPs, equal), and the
physics of the tensors (symmetry-equivalent atoms carry rotated tensors with
the same principal values; B = 8 pi^2 U_eq).

    DIFFUSE_PYTHON_PATH=<build>/python pytest tests/python/test_wp3_adp_molecules.py
"""

from __future__ import annotations

import numpy as np
import pytest

from test_wp3_structures import T

CMM2 = T / "01.Simulation/cmm2_scat.cell"          # water molecules
HEXAGONAL = T / "05.Stack/STACK_TYPES/hexagonal.cell"   # anisotropic Co, P6

ANATASE_ANIS = """\
title Anatase, anisotropic ADPs (test)
spcgr I41/amd
cell 3.785, 3.785, 9.514, 90.0, 90.0, 90.0
scat  Ti, O
adp   0.5, 0.6
occ   1.0, 1.0
anis type:1, values:[0.004, 0.004, 0.006, 0.0, 0.0, 0.0]
anis type:2, values:[0.003, 0.008, 0.005, 0.0, 0.0, 0.0]
atoms
Ti   0.0, 0.0, 0.0,    0.5
O    0.0, 0.0, 0.2081, 0.6
"""
U_TI = [0.004, 0.004, 0.006, 0.0, 0.0, 0.0]
U_O = [0.003, 0.008, 0.005, 0.0, 0.0, 0.0]


@pytest.fixture
def dc(suite):
    import diffuse
    return diffuse


@pytest.fixture
def anatase(dc, tmp_path):
    f = tmp_path / "anatase_anis.cell"
    f.write_text(ANATASE_ANIS)
    return dc.read_structure(f, supercell=(4, 4, 2))


def cartesian_u(crystal, u6):
    """U (CIF convention: U_ij on the reciprocal-length-scaled axes) as a
    Cartesian 3x3 tensor in Å²."""
    a = crystal.lattice_vectors()                  # rows: a, b, c
    n = np.diag(np.linalg.norm(np.linalg.inv(a).T, axis=1))   # |a*|, |b*|, |c*|
    u11, u22, u33, u23, u13, u12 = u6
    u = np.array([[u11, u12, u13], [u12, u22, u23], [u13, u23, u33]])
    return a.T @ n @ u @ n @ a


# ------------------------------------------------------------------ ADPs out
def test_symmetry_rotates_the_tensors(dc, anatase):
    """The four-fold axis of I4_1/amd swaps U11 and U22 of the oxygen site;
    titanium (site -4m2) keeps U11 = U22. Every atom, also beyond the first
    unit cells that DISCUS checks explicitly, has one of these tensors."""
    names, _xyz = anatase.atoms()
    u = anatase.adps()
    assert {tuple(r) for r in u[names == "TI"].round(6)} == {tuple(U_TI)}
    o_swapped = [U_O[1], U_O[0]] + U_O[2:]
    assert {tuple(r) for r in u[names == "O"].round(6)} == {tuple(U_O), tuple(o_swapped)}


def test_principal_values_are_kept(dc, anatase):
    """Rotation keeps the eigenvalues (mean-square displacements)."""
    names, _ = anatase.atoms()
    for name, ref in (("TI", U_TI), ("O", U_O)):
        want = np.linalg.eigvalsh(cartesian_u(anatase, ref))
        for u6 in anatase.adps()[names == name][:20]:
            np.testing.assert_allclose(np.linalg.eigvalsh(cartesian_u(anatase, u6)), want, atol=1e-12)


def test_b_from_the_tensor(dc, anatase):
    """DISCUS sets B = 8 pi^2 U_eq (U_eq: a third of the trace)."""
    names, b, _occ = anatase.atom_types()
    for name, ref in (("TI", U_TI), ("O", U_O)):
        assert b[list(names).index(name)] == pytest.approx(8 * np.pi**2 * np.mean(ref[:3]))


def test_isotropic_atoms_have_the_tensor_of_b(dc):
    cr = dc.read_structure(CMM2, supercell=(2, 1, 1))      # B = 0.1 (O), 0.2 (H), cubic
    names, _ = cr.atoms()
    u = cr.adps()
    for name, b in (("O", 0.1), ("H", 0.2)):
        np.testing.assert_allclose(u[names == name][:, :3], b / (8 * np.pi**2), rtol=1e-9)
        np.testing.assert_allclose(u[names == name][:, 3:], 0.0, atol=1e-15)


# ------------------------------------------------------------------- ADPs in
def test_anisotropic_round_trip(dc, anatase):
    """File -> adps() -> make_crystal(u_aniso per atom) -> adps(): the
    same tensors, and the same Fourier map (which the ADPs change)."""
    names, xyz = anatase.atoms()
    u = anatase.adps()
    kw = dict(ll=(0, 0, 0), lr=(6, 0, 0), ul=(0, 0, 6), npoints=(31, 31))
    ref = dc.fourier(anatase, **kw)
    made = dc.make_crystal(anatase.lattice, names, xyz, u_aniso=u, cells=anatase.supercell,
                           space_group=anatase.space_group)
    np.testing.assert_allclose(made.adps(), u, atol=1e-12)
    np.testing.assert_allclose(made.atom_types()[1], anatase.atom_types()[1][[0, 1, 1]], rtol=1e-9)
    got = dc.fourier(made, **kw)
    np.testing.assert_allclose(got.intensity, ref.intensity, rtol=1e-6, atol=1e-9 * ref.intensity.max())
    iso = dc.make_crystal(anatase.lattice, names, xyz, cells=anatase.supercell,
                          b_iso=dict(zip(*anatase.atom_types()[:2])))
    assert not np.allclose(dc.fourier(iso, **kw).intensity, ref.intensity, rtol=1e-3)


def test_u_aniso_per_name(dc):
    fcc = np.array([[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    u = (0.01, 0.02, 0.03, 0.0, 0.0, 0.0)
    cr = dc.make_crystal((3.6, 3.6, 3.6, 90, 90, 90), ["CU"] * 4, fcc, u_aniso={"cu": u})
    np.testing.assert_allclose(cr.adps(), [u] * 4)
    assert cr.atom_types()[1][0] == pytest.approx(8 * np.pi**2 * 0.02)


def test_hexagonal_tensor_survives(dc):
    """P6 with an anisotropic Co (U12 = U11/2, as the hexagonal metric needs)."""
    cr = dc.read_structure(HEXAGONAL, supercell=(2, 2, 1))
    np.testing.assert_allclose(cr.adps(), [[0.005, 0.005, 0.007, 0.0, 0.0, 0.0025]] * cr.natoms)


# ------------------------------------------------------------------ molecules
def test_molecules_out(dc):
    cr = dc.read_structure(CMM2, supercell=(2, 1, 1))
    types, mols = cr.molecules()
    names, _ = cr.atoms()
    assert len(mols) == 8 and set(types) == {1}
    assert all(sorted(names[m]) == ["H", "H", "O"] for m in mols)
    assert sorted(np.concatenate(mols).tolist()) == list(range(cr.natoms))
    _pos, _types, prop, mole = dc.session()._c.atoms()
    for k, m in enumerate(mols, start=1):
        assert set(mole[m]) == {k}
    assert set(prop) == {3}                    # normal atom + in a molecule


def test_molecules_round_trip(dc):
    """File -> molecules() -> make_crystal(molecules=...) -> molecules(), and
    the same property flags and molecule numbers per atom."""
    cr = dc.read_structure(CMM2, supercell=(2, 2, 1))
    types, mols = cr.molecules()
    names, xyz = cr.atoms()
    _p, _t, prop_ref, mole_ref = dc.session()._c.atoms()
    b = dict(zip(*cr.atom_types()[:2]))
    made = dc.make_crystal(cr.lattice, names, xyz, b_iso=b, cells=cr.supercell,
                           molecules=mols, molecule_types=types)
    types2, mols2 = made.molecules()
    assert types2.tolist() == types.tolist()
    assert [m.tolist() for m in mols2] == [m.tolist() for m in mols]
    _p, _t, prop, mole = dc.session()._c.atoms()
    np.testing.assert_array_equal(prop, prop_ref)
    np.testing.assert_array_equal(mole, mole_ref)


def test_molecules_survive_a_reload(dc):
    names = ["O", "H", "H", "O", "H", "H"]
    xyz = [[0, 0, 0], [0.1, 0, 0], [0, 0.1, 0], [0.5, 0.5, 0.5], [0.6, 0.5, 0.5], [0.5, 0.6, 0.5]]
    w = dc.make_crystal((10, 10, 10, 90, 90, 90), names, xyz, molecules=[[0, 1, 2], [3, 4, 5]],
                        molecule_types=[1, 2])
    dc.read_structure(CMM2)                                  # something else
    types, mols = w.molecules()                              # sends w again
    assert types.tolist() == [1, 2] and [m.tolist() for m in mols] == [[0, 1, 2], [3, 4, 5]]


@pytest.mark.parametrize("kwargs, match", [
    (dict(molecules=[[0, 1], [1, 2]]), "more than one molecule"),
    (dict(molecules=[[0, 7]]), "0 .. natoms-1"),
    (dict(molecules=[[]]), "at least one atom"),
    (dict(molecules=[[0], [1]], molecule_types=[1]), "one type"),
    (dict(molecule_types=[1]), "without molecules"),
    (dict(u_aniso=np.zeros((3, 5))), "shape"),
])
def test_bad_molecules_and_adps(dc, kwargs, match):
    with pytest.raises(ValueError, match=match):
        dc.make_crystal((5, 5, 5, 90, 90, 90), ["O", "H", "H"],
                        [[0, 0, 0], [0.1, 0, 0], [0, 0.1, 0]], **kwargs)


def test_type_first_seen_late_gets_its_tensor(dc):
    """A type that appears only after the atoms DISCUS checks first (5 per
    atom of the unit cell) still gets its own tensor."""
    fcc = np.array([[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    cells = np.array([[i, j, k] for i in range(5) for j in range(5) for k in range(5)])
    xyz = (cells[:, None, :] + fcc[None, :, :]).reshape(-1, 3)
    names = ["CU"] * len(xyz)
    names[-1] = "AU"
    u = np.tile([0.01, 0.01, 0.01, 0, 0, 0], (len(xyz), 1))
    u[-1] = [0.02, 0.03, 0.04, 0, 0, 0]
    cr = dc.make_crystal((3.6, 3.6, 3.6, 90, 90, 90), names, xyz, u_aniso=u, cells=(5, 5, 5))
    np.testing.assert_allclose(cr.adps(), u, atol=1e-12)


@pytest.mark.xfail(strict=True, reason=(
    "DISCUS bug: 'save' writes one atom type per distinct ADP tensor (TI, O, O with "
    "U11 and U22 swapped), but the atom lines only carry the name, so 'read stru' "
    "puts every O on the first O type: the symmetry-rotated tensors are lost "
    "(236 of 384 atoms in anatase 4x4x2). Remove this mark once fixed."))
def test_save_and_read_keep_rotated_tensors(dc, anatase, tmp_path):
    u_ref = anatase.adps()
    dc.run_macro("discus\n  save\n    outfile saved.stru\n    write all\n    run\n  exit\nexit",
                 cwd=tmp_path)
    back = dc.read_structure(tmp_path / "saved.stru")
    np.testing.assert_allclose(back.adps(), u_ref, atol=1e-9)


def test_ase_round_trip_keeps_tensors(dc, anatase):
    pytest.importorskip("ase")
    back = dc.from_ase(anatase.to_ase())
    np.testing.assert_allclose(back.adps(), anatase.adps(), atol=1e-12)
