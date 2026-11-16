"""WP3 safety net: structure reading, space-group expansion, supercells, export/import.

Golden data for the WP3 work package (Python Crystal model and I/O). The
references are independent of DISCUS where possible:

* space-group expansion is compared with gemmi's symmetry operations
  (tests are skipped if gemmi is not installed);
* export/import round trips must give back the same atoms;
* supercells must contain every unit-cell site once per cell.

    DIFFUSE_PYTHON_PATH=<build>/python pytest tests/python/test_wp3_structures.py
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[2]
T = REPO / "tests"

# (file, space group as written in the file, atoms in the unit cell)
STRUCTURES = [
    (T / "02.Calculation_Scattering/cuo.cell", "C2/c", 8),               # monoclinic
    (T / "14.Decoration/Deco_01/CELL/bgamma.cell", "Cmcm", 12),
    (T / "14.Decoration/Deco_01/CELL/boehmite.cell", "Cmc21", 12),
    (T / "08.Nanoparticle_example/COMPLEX/CELL/znse_zincblende.cell", "F-43m", 8),
    (T / "02.Calculation_Scattering/cu2o.cell", "Pn-3m", 6),             # origin choice 1
    (T / "01.Simulation/cmm2_scat.cell", "Cmm2", 12),                     # with molecules
    (T / "01.Simulation/p4mm.cell", "P4mm", 8),
    (T / "03.Pair_Distribution_Function/ni.cell", "Fm-3m", 4),
    (T / "03.Pair_Distribution_Function/water_hex.cell", "P63/mmc", 4),  # hexagonal
]
ANATASE = REPO / "python" / "examples" / "anatase.cif"
CUO = STRUCTURES[0][0]


@pytest.fixture
def dc(suite):
    import diffuse
    return diffuse


# ------------------------------------------------------------------ helpers
def parse_cell(path: Path) -> tuple[str, list[tuple[str, list[float]]]]:
    """Space group and asymmetric unit of a DISCUS .cell file."""
    spcgr, atoms, in_atoms = None, [], False
    for raw in path.read_text().splitlines():
        line = raw.split("#")[0].split("!")[0].strip()
        if not line:
            continue
        words = line.replace(",", " ").split()
        key = words[0].lower()
        if key == "spcgr":
            spcgr = words[1]
        elif key == "atoms":
            in_atoms = True
        elif in_atoms and len(words) >= 4 and key != "molecule":
            try:
                atoms.append((words[0].upper(), [float(v) for v in words[1:4]]))
            except ValueError:
                continue
    return spcgr, atoms


def wrap(d):
    """Map fractional differences into [-0.5, 0.5)."""
    return np.mod(np.asarray(d) + 0.5, 1.0) - 0.5


def element(name: str) -> str:
    """Element symbol from a DISCUS atom name: 'CU2+' -> 'CU'."""
    return re.match(r"[A-Za-z]+", name).group(0).upper()


def assert_same_sites(a, b, *, names=lambda n: n.upper(), tol=1e-4):
    """Two lists of (name, xyz) describe the same set of sites (mod 1)."""
    assert len(a) == len(b), f"{len(a)} sites vs {len(b)}"
    free = list(range(len(b)))
    for name, x in a:
        for k in free:
            m, y = b[k]
            if names(m) == names(name) and np.allclose(wrap(np.subtract(x, y)), 0, atol=tol):
                free.remove(k)
                break
        else:
            raise AssertionError(f"site {name} {np.round(x, 5)} has no partner")


def sites(crystal):
    names, xyz = crystal.atoms()
    return list(zip(names, np.mod(xyz, 1.0)))


def gemmi_expand(spcgr: str, asym):
    gemmi = pytest.importorskip("gemmi")
    group = gemmi.find_spacegroup_by_name(spcgr)
    assert group is not None, f"gemmi does not know {spcgr}"
    out = []
    for name, xyz in asym:
        found = []
        for op in group.operations():
            v = np.mod(np.array(op.apply_to_xyz(xyz)), 1.0)
            if not any(np.allclose(wrap(v - q), 0, atol=1e-4) for q in found):
                found.append(v)
        out += [(name, v) for v in found]
    return out


# ------------------------------------------------------ space-group expansion
@pytest.mark.parametrize("path, spcgr, natoms", STRUCTURES, ids=lambda v: getattr(v, "stem", None))
def test_space_group_expansion_matches_gemmi(dc, path, spcgr, natoms):
    """DISCUS's unit cell = gemmi's symmetry applied to the asymmetric unit."""
    file_spcgr, asym = parse_cell(path)
    assert file_spcgr == spcgr
    cr = dc.read_structure(path)
    assert cr.natoms == natoms
    assert_same_sites(sites(cr), gemmi_expand(spcgr, asym))


def test_cif_expansion_matches_gemmi(dc):
    gemmi = pytest.importorskip("gemmi")
    st = gemmi.read_small_structure(str(ANATASE))
    ref = [(s.type_symbol, np.mod(np.array(s.fract.tolist()), 1.0))
           for s in st.get_all_unit_cell_sites()]
    cr = dc.read_structure(ANATASE)
    assert cr.natoms == 12
    assert_same_sites(sites(cr), ref, names=element)
    cell = st.cell
    assert cr.lattice == pytest.approx((cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma))


HCP_HEAD = "title hcp\nspcgr P63mc\ncell 3.0 3.0 4.8989797 90.0 90.0 120.0\n"


def test_cell_comma_format_evaluates_expressions(dc, tmp_path):
    """Comma format with a column header: atom values go through the
    interpreter, so '1./3.' is one third."""
    f = tmp_path / "hcp.cell"
    f.write_text(HCP_HEAD + "atoms x, y, z, Biso\nco 1./3., 2./3., 0.0, 1.0\n")
    _, xyz = dc.read_structure(f).atoms()
    assert_same_sites([("CO", [1 / 3, 2 / 3, 0.0]), ("CO", [2 / 3, 1 / 3, 0.5])],
                      [("CO", x) for x in np.mod(xyz, 1.0)])


def test_cell_legacy_format_misreads_fractions(dc):
    """Legacy format (bare 'atoms', values separated by blanks): no
    expressions; Fortran list-directed input stops at '/', so '1./3.' is
    read as 1. without a warning. 01.Simulation/calculate.cell uses this
    format and puts Co at (0,0,0) instead of hcp's (1/3,2/3,0). Documents
    the current behaviour; change this test if the reader learns
    expressions or rejects such lines."""
    cr = dc.read_structure(T / "01.Simulation/calculate.cell")
    _, xyz = cr.atoms()
    assert np.allclose(np.mod(xyz[:, :2], 1.0), 0.0)


# ------------------------------------------------------------- supercells
@pytest.mark.parametrize("supercell", [(2, 2, 2), (3, 1, 2), (4, 4, 1)])
def test_supercell_repeats_every_site_once_per_cell(dc, supercell):
    unit = dc.read_structure(CUO)
    unit_sites = sites(unit)
    big = dc.read_structure(CUO, supercell=supercell)
    ncells = int(np.prod(supercell))
    assert big.natoms == ncells * unit.natoms
    names, xyz = big.atoms()
    cells = np.floor(xyz + 1e-6).astype(int)               # which cell each atom is in
    assert len({tuple(c) for c in cells}) == ncells
    for axis, n in enumerate(supercell):
        assert np.ptp(cells[:, axis]) == n - 1             # n consecutive cells per axis
    for name, x in unit_sites:                            # every site once in every cell
        hits = [c for m, y, c in zip(names, xyz, cells)
                if m == name and np.allclose(wrap(y - x), 0, atol=1e-5)]
        assert len({tuple(c) for c in hits}) == ncells and len(hits) == ncells


# --------------------------------------------------------- save / export
def test_save_and_read_back(dc, tmp_path):
    cr = dc.read_structure(T / "01.Simulation/cmm2_scat.cell", supercell=(2, 2, 1))
    dc.run_macro("discus\n  save\n    outfile saved.stru\n    write all\n    run\n  exit\nexit",
                 cwd=tmp_path)
    back = dc.read_structure(tmp_path / "saved.stru")
    assert back.natoms == cr.natoms
    assert back.lattice == pytest.approx(cr.lattice)
    n1, x1 = cr.atoms()
    n2, x2 = back.atoms()
    assert list(n1) == list(n2)
    np.testing.assert_allclose(x1, x2, atol=1e-6)


def export_import(dc, path, fmt, ext, tmp_path):
    src = dc.read_structure(path)
    ref = sites(src)
    dc.run_macro(f"discus\n  export {fmt}, rt.{ext}\nexit", cwd=tmp_path)
    dc.run_macro(f"discus\n  import {fmt}, rt.{ext}\nexit", cwd=tmp_path)
    made = [p for p in tmp_path.iterdir() if p.suffix in (".stru", ".cell")]
    assert len(made) == 1, f"import made {made}"
    return src, ref, dc.read_structure(made[0])


@pytest.mark.parametrize("path", [CUO, ANATASE, T / "01.Simulation/cmm2_scat.cell",
                                  T / "02.Calculation_Scattering/cu2o.cell"],
                         ids=lambda p: p.stem)
def test_cif_export_import_round_trip(dc, path, tmp_path):
    src, ref, back = export_import(dc, path, "cif", "cif", tmp_path)
    assert back.lattice == pytest.approx(src.lattice, abs=1e-4)
    assert_same_sites(ref, sites(back))


def test_exported_cif_read_by_gemmi(dc, tmp_path):
    """An independent reader sees the same unit cell in DISCUS's CIF."""
    gemmi = pytest.importorskip("gemmi")
    cr = dc.read_structure(CUO)
    dc.run_macro("discus\n  export cif, out.cif\nexit", cwd=tmp_path)
    st = gemmi.read_small_structure(str(tmp_path / "out.cif"))
    ref = [(s.type_symbol, np.mod(np.array(s.fract.tolist()), 1.0))
           for s in st.get_all_unit_cell_sites()]
    assert_same_sites(sites(cr), ref, names=element)


@pytest.mark.parametrize("path", [CUO, ANATASE, T / "01.Simulation/cmm2_scat.cell"],
                         ids=lambda p: p.stem)
def test_rmcprofile_export_import_round_trip(dc, path, tmp_path):
    """RMCProfile stores elements without charges ('Cu' for 'CU2+')."""
    _, ref, back = export_import(dc, path, "rmcprofile", "rmc6f", tmp_path)
    assert_same_sites(ref, sites(back), names=element)


def test_vasp_export(dc, tmp_path):
    """VASP 4 layout: lattice vectors, counts per element, direct coordinates."""
    cr = dc.read_structure(CUO)
    dc.run_macro("discus\n  export vasp, out.vasp\nexit", cwd=tmp_path)
    lines = (tmp_path / "out.vasp").read_text().splitlines()
    scale = float(lines[1])
    vectors = scale * np.array([[float(v) for v in lines[i].split()] for i in (2, 3, 4)])
    counts = [int(v) for v in lines[5].split()]
    assert lines[6].strip().lower().startswith("direct")
    xyz = np.array([[float(v) for v in line.split()[:3]] for line in lines[7:7 + sum(counts)]])
    a, b, c = np.linalg.norm(vectors, axis=1)
    alpha = np.degrees(np.arccos(vectors[1] @ vectors[2] / (b * c)))
    beta = np.degrees(np.arccos(vectors[0] @ vectors[2] / (a * c)))
    gamma = np.degrees(np.arccos(vectors[0] @ vectors[1] / (a * b)))
    assert (a, b, c, alpha, beta, gamma) == pytest.approx(cr.lattice, abs=1e-4)
    names, ref = cr.atoms()
    assert counts == [int((names == n).sum()) for n in dict.fromkeys(names)]
    np.testing.assert_allclose(xyz, ref, atol=1e-6)


@pytest.mark.xfail(strict=True, reason=(
    "export shelx writes every atom of the unit cell but keeps the space-group "
    "symmetry (LATT/SYMM) and sets occupancies to 0.5; re-importing CuO gives "
    "32 atoms instead of 8. Remove this mark once the export is fixed."))
def test_shelx_export_import_round_trip(dc, tmp_path):
    _, ref, back = export_import(dc, CUO, "shelx", "ins", tmp_path)
    assert_same_sites(ref, sites(back), names=element)


def test_shelx_export_refuses_origin_choice_2(dc, tmp_path):
    """SHELX needs the centre of symmetry at the origin; DISCUS says so."""
    dc.read_structure(ANATASE)
    with pytest.raises(dc.MacroError, match="1bar not at origin"):
        dc.run_macro("discus\n  export shelx, out.ins\nexit", cwd=tmp_path)
