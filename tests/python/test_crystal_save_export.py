"""Crystal.save() and Crystal.export(), each checked by reading back what it
wrote: with DISCUS (read_structure, import) and, where possible, with an
independent reader (gemmi for CIF) or by parsing the file.

    DIFFUSE_PYTHON_PATH=<build>/python pytest tests/python/test_crystal_save_export.py
"""

from __future__ import annotations

import numpy as np
import pytest

from test_wp3_structures import ANATASE, T, assert_same_sites, element, sites

CMM2 = T / "01.Simulation/cmm2_scat.cell"            # molecules
CUO = T / "02.Calculation_Scattering/cuo.cell"        # monoclinic, two types
HEXAGONAL = T / "05.Stack/STACK_TYPES/hexagonal.cell"   # anisotropic tensor
NICKEL = T / "03.Pair_Distribution_Function/ni.cell"


@pytest.fixture
def dc(suite):
    import diffuse
    return diffuse


def snapshot(dc, crystal):
    """Everything about the crystal DISCUS holds, through the C interface."""
    names, xyz = crystal.atoms()
    tnames, b, occ = crystal.atom_types()
    mtypes, mols = crystal.molecules()
    _pos, types, prop, mole = dc.session()._c.atoms()
    lattice, spcgr, cells, per_cell = dc.session()._c.cell()
    return dict(names=names.tolist(), xyz=xyz, types=types.tolist(), tnames=tnames.tolist(),
                b=b, occ=occ, adps=crystal.adps(), prop=prop.tolist(), mole=mole.tolist(),
                mtypes=mtypes.tolist(), mols=[m.tolist() for m in mols],
                lattice=lattice, spcgr=spcgr, cells=cells, per_cell=per_cell)


def assert_same(a, b):
    for key in a:
        if isinstance(a[key], np.ndarray):
            np.testing.assert_allclose(b[key], a[key], atol=1e-6, err_msg=key)
        elif key == "lattice":
            assert b[key] == pytest.approx(a[key], abs=1e-6), key
        else:
            assert b[key] == a[key], key


# ---------------------------------------------------------------------- save
@pytest.mark.parametrize("suffix", [".stru", ".h5"])
@pytest.mark.parametrize("path, supercell", [(CMM2, (2, 2, 1)), (CUO, (2, 2, 2)),
                                             (HEXAGONAL, (3, 3, 2))],
                         ids=["molecules", "monoclinic", "anisotropic"])
def test_save_and_read_back(dc, tmp_path, path, supercell, suffix):
    """Everything comes back: atoms, types, B, occupancy, tensors, property
    flags, molecules, cell, space group, number of unit cells."""
    cr = dc.read_structure(path, supercell=supercell)
    ref = snapshot(dc, cr)
    written = cr.save(tmp_path / f"saved{suffix}")
    assert written == tmp_path / f"saved{suffix}" and written.is_file()
    back = dc.read_structure(written)
    assert back.natoms == cr.natoms
    assert_same(ref, snapshot(dc, back))


def made_crystal(dc, cells, u_aniso=None):
    """Cu with Au on every fourth site, made in Python, with two molecules."""
    molecules = [[0, 1], [4, 5, 6]] if np.prod(cells) > 1 else [[0, 1], [2, 3]]
    fcc = np.array([[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    grid = np.array([[i, j, k] for i in range(cells[0]) for j in range(cells[1])
                     for k in range(cells[2])])
    xyz = (grid[:, None, :] + fcc[None]).reshape(-1, 3)
    names = np.where(np.arange(len(xyz)) % 4 == 0, "AU", "CU")
    return dc.make_crystal((3.75, 3.75, 3.75, 90, 90, 90), names, xyz, space_group="P1",
                           b_iso={"AU": 0.6, "CU": 0.8}, u_aniso=u_aniso, cells=cells,
                           molecules=molecules, molecule_types=[1, 2])


ANISO = {"AU": (0.01, 0.01, 0.01, 0, 0, 0), "CU": (0.008, 0.03, 0.01, 0, 0, 0)}


def test_save_a_crystal_made_in_python(dc, tmp_path):
    """make_crystal (B per name, molecules, a 3 x 2 x 2 supercell) -> save ->
    read: everything, and the same Fourier map. (The per-atom tensors that
    adps() reports differ after 'read stru' beyond the first 20 atoms, see
    test_save_anisotropic_supercell; isotropic calculations use B per type
    and are not affected.)"""
    cr = made_crystal(dc, (3, 2, 2))
    ref = snapshot(dc, cr)
    kw = dict(ll=(0, 0, 0), lr=(6, 0, 0), ul=(0, 6, 0), npoints=(31, 31))
    f_ref = dc.fourier(cr, **kw).intensity
    back = dc.read_structure(cr.save(tmp_path / "made.stru"))
    got = snapshot(dc, back)
    ref.pop("adps"), got.pop("adps")
    assert_same(ref, got)
    np.testing.assert_allclose(dc.fourier(back, **kw).intensity, f_ref, rtol=1e-9)


def test_save_anisotropic_unit_cell(dc, tmp_path):
    """Anisotropic tensors per type, one unit cell: they read back."""
    cr = made_crystal(dc, (1, 1, 1), u_aniso=ANISO)
    ref = snapshot(dc, cr)
    back = dc.read_structure(cr.save(tmp_path / "made.stru"))
    assert_same(ref, snapshot(dc, back))


@pytest.mark.xfail(strict=True, reason=(
    "DISCUS bug: 'read stru' sets the displacement tensors of the first "
    "5 x (atoms per unit cell) atoms only (prep_anis); every later atom keeps tensor 1 "
    "(isotropic crystals only report it wrong in adps(); their calculations use B). "
    "For this 3x2x2 Cu/Au crystal 21 of 48 atoms get the wrong tensor and a Fourier map "
    "changes by up to 3 %. dc.make_crystal works around it; reading a .cell file is not "
    "affected. Remove this mark once fixed."))
def test_save_anisotropic_supercell(dc, tmp_path):
    cr = made_crystal(dc, (3, 2, 2), u_aniso=ANISO)
    u = cr.adps()
    back = dc.read_structure(cr.save(tmp_path / "made.stru"))
    np.testing.assert_allclose(back.adps(), u, atol=1e-9)


def test_save_reloads_an_older_crystal(dc, tmp_path):
    """Saving a crystal that DISCUS no longer holds loads it first."""
    nickel = dc.read_structure(NICKEL, supercell=2)
    dc.read_structure(CUO)
    back = dc.read_structure(nickel.save(tmp_path / "ni.stru"))
    assert back.natoms == 32 and set(back.atoms()[0]) == {"NI"}


def test_save_replaces_an_old_file(dc, tmp_path):
    out = tmp_path / "x.stru"
    out.write_text("old")
    dc.read_structure(NICKEL).save(out)
    assert out.read_text().startswith("title")


@pytest.mark.parametrize("name, error", [("a,b.stru", ValueError), ("no/such/dir/x.stru", FileNotFoundError)])
def test_save_bad_names(dc, tmp_path, name, error):
    with pytest.raises(error):
        dc.read_structure(NICKEL).save(tmp_path / name)


def test_save_refuses_export_formats(dc, tmp_path):
    with pytest.raises(ValueError, match="export"):
        dc.read_structure(NICKEL).save(tmp_path / "x.cif", format="cif")


# -------------------------------------------------------------------- export
def supercell_frame(crystal):
    """Fractional coordinates in units of the whole crystal (what export
    writes as one cell)."""
    names, xyz = crystal.atoms()
    return names, xyz / np.array(crystal.supercell, dtype=float)


def test_export_cif_read_by_gemmi(dc, tmp_path):
    """P1, the whole crystal as one cell, every atom, U_iso = B / 8 pi^2."""
    gemmi = pytest.importorskip("gemmi")
    cr = dc.read_structure(CUO, supercell=(2, 2, 1))
    st = gemmi.read_small_structure(str(cr.export(tmp_path / "cuo.cif")))
    a, b, c, al, be, ga = cr.lattice
    assert (st.cell.a, st.cell.b, st.cell.c) == pytest.approx((2 * a, 2 * b, c), abs=1e-4)
    assert (st.cell.alpha, st.cell.beta, st.cell.gamma) == pytest.approx((al, be, ga), abs=1e-4)
    assert st.spacegroup_hm.replace(" ", "") == "P1"
    names, frac = supercell_frame(cr)
    got = [(s.type_symbol, np.array(s.fract.tolist())) for s in st.sites]
    assert len(got) == cr.natoms
    assert_same_sites(list(zip(names, frac)), got, names=element)
    tnames, b_iso, _ = cr.atom_types()
    u = {element(n): bb / (8 * np.pi**2) for n, bb in zip(tnames, b_iso)}
    for s in st.sites:
        assert s.u_iso == pytest.approx(u[element(s.type_symbol)], abs=1e-6)


def test_export_cif_read_back_by_discus(dc, tmp_path):
    cr = dc.read_structure(CUO, supercell=(2, 2, 1))
    names, frac = supercell_frame(cr)
    back = dc.read_structure(cr.export(tmp_path / "cuo.cif"))
    assert back.natoms == cr.natoms
    assert back.lattice[:3] == pytest.approx((2 * cr.lattice[0], 2 * cr.lattice[1], cr.lattice[2]),
                                             abs=1e-4)
    assert_same_sites(list(zip(names, np.mod(frac, 1))), sites(back))


@pytest.mark.xfail(strict=True, reason=(
    "DISCUS bug: 'export cif' with 'spcgr:original' writes every atom of the unit cell "
    "together with the space group's operations, so readers (DISCUS and gemmi alike) "
    "generate 80 atoms for anatase instead of 12. Remove this mark once fixed."))
def test_export_cif_original_space_group(dc, tmp_path):
    """For one unit cell, the original space group reads back to the same
    crystal (DISCUS expands it again)."""
    cr = dc.read_structure(ANATASE)
    out = cr.export(tmp_path / "anatase.cif", space_group="original")
    assert "I 41/a m d" in out.read_text() or "I41/amd" in out.read_text().replace(" ", "")
    back = dc.read_structure(out)
    assert back.natoms == cr.natoms
    assert_same_sites(sites(cr), sites(back))


def test_export_rmcprofile_read_back(dc, tmp_path):
    """DISCUS's own import reads the file back: the same atoms in the same
    supercell (RMCProfile keeps elements, not charges)."""
    cr = dc.read_structure(ANATASE, supercell=(2, 2, 1))
    out = cr.export(tmp_path / "anatase.rmc6f", version=6)
    text = out.read_text()
    assert f"Number of atoms:{cr.natoms:>26d}" in text or f"Number of atoms: " in text
    dc.run_macro(f"discus\n  import rmcprofile, {out.name}\nexit", cwd=tmp_path)
    made = [p for p in tmp_path.iterdir() if p.suffix == ".stru"]
    assert len(made) == 1
    back = dc.read_structure(made[0])
    back_cells = dc.session()._c.cell()[2]          # while DISCUS holds it
    assert back.natoms == cr.natoms
    assert sorted(element(n) for n in back.atoms()[0]) == sorted(element(n) for n in cr.atoms()[0])
    # The export moves the crystal to start at 0 (a shift by whole unit
    # cells: the lower corner of the crystal) and the import groups the
    # atoms by element.
    names, xyz = cr.atoms()
    bnames, bxyz = back.atoms()
    bxyz = bxyz * np.array(cr.supercell) / np.array(back_cells)
    box = np.array(cr.supercell, dtype=float)
    moved = xyz - np.floor(xyz.min(axis=0) + 1e-6)
    assert_same_sites(list(zip(names, np.mod(moved, box) / box)),
                      list(zip(bnames, np.mod(bxyz, box) / box)), names=element)


def test_export_vasp_parsed(dc, tmp_path):
    """POSCAR: the whole crystal's lattice vectors, counts per element and
    direct coordinates of every atom."""
    cr = dc.read_structure(CUO, supercell=(1, 2, 2))
    out = cr.export(tmp_path / "POSCAR")
    lines = out.read_text().splitlines()
    vectors = float(lines[1]) * np.array([[float(v) for v in lines[i].split()] for i in (2, 3, 4)])
    counts = [int(v) for v in lines[5].split()]
    assert lines[6].strip().lower().startswith("direct")
    xyz = np.array([[float(v) for v in ln.split()[:3]] for ln in lines[7:7 + sum(counts)]])
    np.testing.assert_allclose(vectors, cr.lattice_vectors() * np.array([[1], [2], [2]]), atol=1e-6)
    names, frac = supercell_frame(cr)
    assert counts == [int((names == n).sum()) for n in dict.fromkeys(names)]
    assert_same_sites(list(zip(names, np.mod(frac, 1))), list(zip(np.repeat(list(dict.fromkeys(names)), counts), np.mod(xyz, 1))))


def test_export_scatty_parsed(dc, tmp_path):
    """Scatty: the unit cell, the box of cells, the sites, and per atom its
    site and cell; site + cell gives back every position."""
    cr = dc.read_structure(CUO, supercell=(2, 2, 1))
    out = cr.export(tmp_path / "cuo_atoms_01.txt", format="scatty")
    cell, box, site_xyz, atoms = None, None, [], []
    for line in out.read_text().splitlines():
        w = line.split()
        if not w:
            continue
        if w[0] == "CELL":
            cell = [float(v) for v in w[1:7]]
        elif w[0] == "BOX":
            box = [int(v) for v in w[1:4]]
        elif w[0] == "SITE":
            site_xyz.append([float(v) for v in w[1:4]])
        elif w[0] == "ATOM":
            atoms.append((w[-1], int(w[1]), [int(v) for v in w[2:5]]))
    assert cell == pytest.approx(list(cr.lattice), abs=1e-4) and box == [2, 2, 1]
    assert len(atoms) == cr.natoms
    pos = [(name, np.mod(np.array(site_xyz[s - 1]) + ijk, box) / box) for name, s, ijk in atoms]
    names, frac = supercell_frame(cr)
    assert_same_sites(list(zip(names, np.mod(frac, 1))), pos, names=element)


def test_export_shelx_writes_a_file(dc, tmp_path):
    """The SHELX file is written with the cell and every atom; reading it back
    hits a known DISCUS bug (test_wp3_structures.test_shelx_export_import_round_trip)."""
    cr = dc.read_structure(CUO)
    text = cr.export(tmp_path / "cuo.ins", cycles=5, wavelength="MOA1").read_text()
    assert "CELL" in text and "L.S." in text


@pytest.mark.parametrize("name, kwargs, error, match", [
    ("x.xyz", {}, ValueError, "format"),
    ("x.txt", {"format": "cif"}, ValueError, "suffix"),
    ("x.rmc6f", {"version": 5}, ValueError, "version"),
    ("x_atoms_01.txt", {"format": "scatty", "site": "every"}, ValueError, "site"),
    ("x.cif", {"cycles": 3}, TypeError, "SHELX"),
    ("x.ins", {"colour": "red"}, TypeError, "unknown"),
])
def test_export_bad_arguments(dc, tmp_path, name, kwargs, error, match):
    with pytest.raises(error, match=match):
        dc.read_structure(NICKEL).export(tmp_path / name, **kwargs)
