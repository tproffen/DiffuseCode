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
    assert capi["version"] == 1
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


def test_reset(capi):
    assert capi["atoms_before"] == 32           # fcc, 2x2x2 cells
    assert capi["atoms_after_reset"] == 0
