"""The C interface (libdiffuse_core, python/diffuse_core.h), called with ctypes.

Every function of the header is called once, straight through ctypes and
without the diffuse package, in a separate process (one DISCUS per
process). Option B in migration/WP2_Binding_Decision.md.

    DIFFUSE_PYTHON_PATH=<build>/python pytest tests/python/test_capi.py
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[2]
NAME = {"darwin": "libdiffuse_core.dylib"}.get(sys.platform, "libdiffuse_core.so")


def find_library() -> Path | None:
    """The library of the build in use: $DIFFUSE_PYTHON_PATH, else the
    installed diffuse package."""
    folders = []
    if os.environ.get("DIFFUSE_PYTHON_PATH"):
        folders.append(Path(os.environ["DIFFUSE_PYTHON_PATH"]))
    else:
        import importlib.util
        spec = importlib.util.find_spec("diffuse")
        if spec and spec.origin:
            folders.append(Path(spec.origin).parent)
    for folder in folders:
        if (folder / NAME).exists():
            return folder / NAME
    return None


SCRIPT = r'''
import ctypes, json, os, sys
import numpy as np
lib_path, work = sys.argv[1], sys.argv[2]
os.chdir(work)
real_stdout = os.dup(1)
os.dup2(os.open(os.devnull, os.O_WRONLY), 1)      # DISCUS output: not needed here
lib = ctypes.CDLL(lib_path)
i32, i64, f64, ch = ctypes.c_int, ctypes.c_int64, ctypes.c_double, ctypes.c_char_p
P = ctypes.POINTER
f32a = np.ctypeslib.ndpointer(np.float32, flags="C_CONTIGUOUS")
i32a = np.ctypeslib.ndpointer(np.int32, flags="C_CONTIGUOUS")
i64a = np.ctypeslib.ndpointer(np.int64, flags="C_CONTIGUOUS")
f64a = np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS")
sig = {
    "dc_api_version": ([], i32), "dc_init": ([], None), "dc_reset": ([], None),
    "dc_flush_output": ([], None), "dc_run": ([ch, i64], i64),
    "dc_value": ([ch, i64], f64), "dc_set_value": ([ch, i64, f64], None),
    "dc_clear_error": ([], None),
    "dc_last_error": ([P(i32), P(i32), ch, ch, ch], None),
    "dc_seed": ([i32], None), "dc_random_size": ([], i64),
    "dc_get_random_state": ([i32a, i64], i32), "dc_set_random_state": ([i32a, i64], i32),
    "dc_data_sets": ([], i32), "dc_data_length": ([i32], i64),
    "dc_get_data": ([i32, f32a, f32a, i64], i32),
    "dc_get_data_2d": ([i32, f32a, f32a, f32a, i64, P(i64), P(i64)], i32),
    "dc_crystal_size": ([P(i64), P(i32)], None),
    "dc_get_cell": ([f64a, ch, i64a, P(i64)], None),
    "dc_get_atoms": ([i64, f64a, i32a, i32a, i32a], i32),
    "dc_get_atom_types": ([i32, ch, f64a, f64a], i32),
    "dc_set_crystal": ([f64a, ch, i64, i64a, i64, i32, ch, f64a, f64a, i32, f64a,
                        i64, f64a, i32a, i32a], i32),
    "dc_get_type_anis": ([i32, f64a], i32), "dc_adp_count": ([], i64),
    "dc_get_adps": ([i64, f64a, f64a], i32), "dc_get_atom_adps": ([i64, i32a], i32),
    "dc_molecule_size": ([P(i64), P(i32), P(i64)], None),
    "dc_get_molecules": ([i64, i64, i32a, i32a, i32a, i32a], i32),
    "dc_set_molecules": ([i64, i64, i32a, i32a, i32a], i32),
}
for name, (args, res) in sig.items():
    f = getattr(lib, name); f.argtypes = args; f.restype = res

def run(text):
    b = text.encode(); return lib.dc_run(b, len(b))
def value(expr):
    b = expr.encode(); return lib.dc_value(b, len(b))
def last_error():
    num, typ = i32(), i32()
    cls, text, msgs = (ctypes.create_string_buffer(n) for n in (4, 45, 560))
    lib.dc_last_error(ctypes.byref(num), ctypes.byref(typ), cls, text, msgs)
    return num.value, cls.raw.decode().strip(), text.raw.decode().strip(), msgs.raw.decode()

r = {}
r["version"] = lib.dc_api_version()
lib.dc_init()
r["run"] = run("r[1] = 6*7")
r["value"] = value("r[1]")
r["third"] = value("1.0/3.0")
run("variable real, myval")
b = b"myval"; lib.dc_set_value(b, len(b), 2.5); r["set_value"] = value("myval")
lib.dc_clear_error(); b = b"r[2]"; lib.dc_set_value(b, len(b), 2.5); r["set_r"] = last_error()[:3]
lib.dc_clear_error(); run("discus\n  no_such_command\nexit"); r["error"] = last_error()[:3]
lib.dc_clear_error(); r["bad_value"] = value("no_such_variable[3]"); r["value_error"] = last_error()[:2]
lib.dc_clear_error(); run("discus\n  read\n    stru /no/such/file.stru\nexit")
r["details"] = "/no/such/file.stru" in last_error()[3]
r["long_line"] = run("echo \"" + "x" * 1100 + "\"")
lib.dc_seed(7); n = lib.dc_random_size(); state = np.zeros(n, np.int32)
r["get_state"] = lib.dc_get_random_state(state, n)
a = [value("ran(0)") for _ in range(3)]
r["set_state"] = lib.dc_set_random_state(state, n)
r["repeat"] = a == [value("ran(0)") for _ in range(3)]
r["state_bad_n"] = lib.dc_get_random_state(np.zeros(n + 1, np.int32), n + 1)
run("kuplot\n  rese\n  func r[0]**2, 0.0, 3.0, 0.5\nexit")
r["sets"] = lib.dc_data_sets(); n1 = lib.dc_data_length(1)
x = np.zeros(n1, np.float32); y = np.zeros(n1, np.float32)
r["get_data"] = lib.dc_get_data(1, x, y, n1); r["xy"] = [x.tolist(), y.tolist()]
r["data_bad_n"] = lib.dc_get_data(1, x, y, n1 - 1); r["no_set"] = lib.dc_data_length(99)
nx, ny = i64(), i64(); z = np.zeros(n1, np.float32)
r["not_a_grid"] = lib.dc_get_data_2d(1, x, y, z, n1, ctypes.byref(nx), ctypes.byref(ny))
cell = os.path.join(sys.argv[3], "python", "examples", "nickel.cell")
run(f"""discus
  read
    cell {cell}, 2, 2, 2
  fourier
    neutron
    ll 0,0,0
    lr 2,0,0
    ul 0,1,0
    na 5
    no 3
    abs h
    ord k
    run
  exit
  output
    outf grid.inte
    value inte
    form stan
    run
  exit
exit
kuplot
  rese
  load ni, grid.inte
exit""")
n2 = lib.dc_data_length(1); gx, gy, gz = (np.zeros(n2, np.float32) for _ in range(3))
r["grid"] = [lib.dc_get_data_2d(1, gx, gy, gz, n2, ctypes.byref(nx), ctypes.byref(ny)), nx.value, ny.value, n2]
r["grid_x"] = sorted(set(np.round(gx, 4).tolist())); r["grid_y"] = sorted(set(np.round(gy, 4).tolist()))
run("discus\n  read\n    cell " + cell + ", 2, 3, 1\nexit")
na, nt = i64(), i32(); lib.dc_crystal_size(ctypes.byref(na), ctypes.byref(nt))
r["size"] = [na.value, nt.value]
lat = np.zeros(6); sg = ctypes.create_string_buffer(16); icc = np.zeros(3, np.int64); nc = i64()
lib.dc_get_cell(lat, sg, icc, ctypes.byref(nc))
r["cell"] = [lat.tolist(), sg.raw.decode().strip(), icc.tolist(), nc.value]
pos = np.zeros((na.value, 3)); ty, pr, mo = (np.zeros(na.value, np.int32) for _ in range(3))
r["get_atoms"] = lib.dc_get_atoms(na.value, pos, ty, pr, mo)
r["pos_range"] = [pos.min(axis=0).tolist(), pos.max(axis=0).tolist()]
r["types"] = sorted(set(ty.tolist()))
r["atoms_bad_n"] = lib.dc_get_atoms(na.value + 1, np.zeros((na.value + 1, 3)),
                                    *(np.zeros(na.value + 1, np.int32) for _ in range(3)))
names = ctypes.create_string_buffer(4 * (nt.value + 1)); b = np.zeros(nt.value + 1); occ = np.zeros(nt.value + 1)
r["get_types"] = lib.dc_get_atom_types(nt.value, names, b, occ)
r["type_table"] = [names.raw.decode(), b.tolist(), occ.tolist()]
r["types_bad_n"] = lib.dc_get_atom_types(nt.value + 1, ctypes.create_string_buffer(4 * (nt.value + 2)),
                                         np.zeros(nt.value + 2), np.zeros(nt.value + 2))
# Back in: a 2 x 1 x 1 NaCl-like crystal from arrays, then read out again
lat = np.array([5.0, 5.0, 5.0, 90.0, 90.0, 90.0]); sg = b"F m -3 m"
pos_in = np.array([[0, 0, 0], [0.5, 0.5, 0.5], [1, 0, 0], [1.5, 0.5, 0.5]], float)
ty_in = np.array([1, 2, 1, 2], np.int32); pr_in = np.ones(4, np.int32)
anis_in = np.array([[0.01, 0.01, 0.01, 0, 0, 0], [0.02, 0.02, 0.02, 0, 0, 0]], float)
r["set"] = lib.dc_set_crystal(lat, sg, len(sg), np.array([2, 1, 1], np.int64), 2, 2, b"NA  CL  ",
                              np.array([0.5, 0.7]), np.array([1.0, 1.0]), 1, anis_in.ravel(),
                              4, pos_in, ty_in, pr_in)
ta = np.zeros((3, 6)); r["type_anis"] = [lib.dc_get_type_anis(2, ta), ta.tolist()]
nan = lib.dc_adp_count(); uij = np.zeros((nan, 6)); prin = np.zeros((nan, 3, 4))
r["adps"] = [nan, lib.dc_get_adps(nan, uij, prin), uij.tolist(), prin[:, :, 3].tolist()]
idx = np.zeros(4, np.int32); r["atom_adps"] = [lib.dc_get_atom_adps(4, idx), idx.tolist()]
r["adps_bad_n"] = lib.dc_get_adps(nan + 1, np.zeros((nan + 1, 6)), np.zeros((nan + 1, 3, 4)))
r["set_mol"] = lib.dc_set_molecules(2, 4, np.array([1, 2], np.int32), np.array([2, 2], np.int32),
                                    np.array([1, 2, 3, 4], np.int32))
nm, nmt, nma = i64(), i32(), i64()
lib.dc_molecule_size(ctypes.byref(nm), ctypes.byref(nmt), ctypes.byref(nma))
mt, mc, ml = (np.zeros(2, np.int32) for _ in range(3)); ma = np.zeros(4, np.int32)
r["mol"] = [nm.value, nmt.value, nma.value, lib.dc_get_molecules(2, 4, mt, mc, ml, ma),
            mt.tolist(), mc.tolist(), ml.tolist(), ma.tolist()]
lib.dc_clear_error()
r["set_mol_bad"] = lib.dc_set_molecules(1, 1, np.array([1], np.int32), np.array([1], np.int32),
                                        np.array([9], np.int32))
lib.dc_crystal_size(ctypes.byref(na), ctypes.byref(nt)); r["set_size"] = [na.value, nt.value]
pos = np.zeros((4, 3)); ty, pr, mo = (np.zeros(4, np.int32) for _ in range(3))
lib.dc_get_atoms(4, pos, ty, pr, mo); r["set_atoms"] = [pos.tolist(), ty.tolist()]
lib.dc_get_cell(lat, (sgb := ctypes.create_string_buffer(16)), icc, ctypes.byref(nc))
r["set_cell"] = [sgb.raw.decode().strip(), icc.tolist(), nc.value]
lib.dc_clear_error()
r["set_bad_type"] = lib.dc_set_crystal(lat, sg, len(sg), np.array([2, 1, 1], np.int64), 2, 2,
                                       b"NA  CL  ", np.zeros(2), np.ones(2), 0, np.zeros(12),
                                       4, pos_in,
                                       np.array([1, 3, 1, 2], np.int32), pr_in)
r["set_bad_error"] = last_error()[:1]
lib.dc_clear_error()
r["set_bad_sg"] = lib.dc_set_crystal(lat, b"Q 99", 4, np.array([2, 1, 1], np.int64), 2, 2,
                                     b"NA  CL  ", np.zeros(2), np.ones(2), 0, np.zeros(12),
                                     4, pos_in, ty_in, pr_in)
r["sg_error"] = last_error()[3]
run("discus\n  read\n    cell " + cell + ", 2, 2, 2\nexit")
run("discus\n  r[3] = n[1]\nexit"); r["atoms_before"] = value("r[3]")
lib.dc_reset(); run("discus\n  r[3] = n[1]\nexit"); r["atoms_after_reset"] = value("r[3]")
lib.dc_flush_output()
os.write(real_stdout, json.dumps(r).encode())
'''


@pytest.fixture(scope="module")
def capi(tmp_path_factory):
    lib = find_library()
    if lib is None:
        if os.environ.get("DIFFUSE_REQUIRE_PYTHON"):
            pytest.fail(f"{NAME} not found next to the Python module")
        pytest.skip(f"{NAME} not found (build with -DDIFFUSE_PYTHON=ON)")
    work = tmp_path_factory.mktemp("capi")
    script = work / "capi_script.py"
    script.write_text(SCRIPT)
    proc = subprocess.run([sys.executable, str(script), str(lib), str(work), str(REPO)],
                          capture_output=True, text=True, timeout=300)
    assert proc.returncode == 0, proc.stdout[-3000:] + proc.stderr[-3000:]
    return json.loads(proc.stdout)


def test_version_and_commands(capi):
    assert capi["version"] == 4
    assert capi["run"] == 0
    assert capi["value"] == 42.0


def test_value_precision(capi):
    """dc_value returns a double: about 12 significant digits, as many as
    DISCUS's expression evaluator keeps (the old f2py module gave about 7)."""
    assert abs(capi["third"] - 1.0 / 3.0) < 1e-11


def test_set_value(capi):
    """dc_set_value sets user variables ("variable real, name"); r[] and
    other built-in arrays are set with commands (dc_run "r[2] = ...")."""
    assert capi["set_value"] == 2.5
    assert capi["set_r"][0] == -24              # "Variable is not defined"


def test_errors(capi):
    assert capi["error"] == [-8, "COMM", "Unknown command"]
    assert capi["bad_value"] == 0.0 and capi["value_error"] == [-2, "FORT"]
    assert capi["details"]


def test_long_line_is_refused(capi):
    assert capi["long_line"] > 1024


def test_random_state(capi):
    assert capi["get_state"] == 0 and capi["set_state"] == 0
    assert capi["repeat"]
    assert capi["state_bad_n"] == -1


def test_kuplot_data(capi):
    assert capi["sets"] == 1 and capi["get_data"] == 0
    x, y = capi["xy"]
    assert x == [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
    assert y == [v * v for v in x]
    assert capi["data_bad_n"] == -1 and capi["no_set"] == -1
    assert capi["not_a_grid"] == -1


def test_grid_data(capi):
    ok, nx, ny, n = capi["grid"]
    assert (ok, nx, ny, n) == (0, 5, 3, 15)
    assert capi["grid_x"] == [0.0, 0.5, 1.0, 1.5, 2.0]
    assert capi["grid_y"] == [0.0, 0.5, 1.0]


def test_crystal_arrays(capi):
    """nickel.cell (F m -3 m, 4 atoms per cell) as 2 x 3 x 1 unit cells."""
    assert capi["size"] == [24, 1]
    lattice, spcgr, icc, ncatoms = capi["cell"]
    assert lattice == [3.5234, 3.5234, 3.5234, 90.0, 90.0, 90.0]
    assert spcgr == "F m -3 m" and icc == [2, 3, 1] and ncatoms == 4
    assert capi["get_atoms"] == 0 and capi["types"] == [1]
    # x spans 2 cells, y 3, z 1: the order (n, 3) of the array is right
    lo, hi = capi["pos_range"]
    assert [h - l for l, h in zip(lo, hi)] == [1.5, 2.5, 0.5]
    assert capi["atoms_bad_n"] == -1
    names, b, occ = capi["type_table"]
    assert capi["get_types"] == 0 and names == "VOIDNI  "
    assert b == [0.0, 0.2] and occ == [1.0, 1.0]
    assert capi["types_bad_n"] == -1


def test_set_crystal(capi):
    assert capi["set"] == 0 and capi["set_size"] == [4, 2]
    assert capi["set_atoms"] == [[[0, 0, 0], [0.5, 0.5, 0.5], [1, 0, 0], [1.5, 0.5, 0.5]],
                                 [1, 2, 1, 2]]
    assert capi["set_cell"] == ["F m -3 m", [2, 1, 1], 2]
    assert capi["set_bad_type"] == -1 and capi["set_bad_error"][0] < 0
    assert capi["set_bad_sg"] == -1 and "Q 99" in capi["sg_error"]


def test_adps_and_molecules(capi):
    """Two isotropic types given as U (one value): two tensors U*1, each atom
    pointing to its type's; molecules in and out."""
    ok, table = capi["type_anis"]
    assert ok == 0 and table[1][0] == 0.01 and table[2][0] == 0.02
    nanis, ok, uij, msd = capi["adps"]
    assert nanis == 2 and ok == 0
    assert [u[:3] for u in uij] == [[0.01] * 3, [0.02] * 3]
    assert np.allclose([sorted(m) for m in msd], [[0.01] * 3, [0.02] * 3])  # <u^2> per axis
    assert capi["atom_adps"] == [0, [1, 2, 1, 2]]
    assert capi["adps_bad_n"] == -1
    assert capi["set_mol"] == 0
    assert capi["mol"] == [2, 2, 4, 0, [1, 2], [0, 0], [2, 2], [1, 2, 3, 4]]
    assert capi["set_mol_bad"] == -1


def test_reset(capi):
    assert capi["atoms_before"] == 32           # fcc, 2x2x2 cells
    assert capi["atoms_after_reset"] == 0
