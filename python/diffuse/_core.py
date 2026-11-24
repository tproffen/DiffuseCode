"""The compiled DISCUS core: the C interface ``libdiffuse_core``
(``python/diffuse_core.h``), loaded with ctypes
(migration/WP2_Binding_Decision.md, option B).

Every function is declared in ``_SIGNATURES`` below, in one place; a wrong
declaration shows up at once in the tests. ``Session`` (``_session.py``)
uses nothing but ``CtypesCore``. Strings go in as ``str``; arrays come back
as NumPy arrays. The library needs no particular Python version.
"""

from __future__ import annotations

import ctypes
import os
import sys
from pathlib import Path

import numpy as np

LIBRARY = {"darwin": "libdiffuse_core.dylib", "win32": "diffuse_core.dll"}.get(
    sys.platform, "libdiffuse_core.so")
API_VERSION = 1              # DC_API_VERSION in diffuse_core.h
ERR_LENGTHS = (4, 45, 560)   # DC_ERR_CLASS_LEN, DC_ERR_TEXT_LEN, DC_ERR_MSGS_LEN


def _folders() -> list[tuple[Path, bool]]:
    """Where to look, in order: (folder, is the installed package).

    $DIFFUSE_PYTHON_PATH (a build folder's python/) wins, for development
    builds; then the copy that "pip install" put inside this package.
    """
    folders = []
    if os.environ.get("DIFFUSE_PYTHON_PATH"):
        folders.append((Path(os.environ["DIFFUSE_PYTHON_PATH"]), False))
    folders.append((Path(__file__).resolve().parent, True))
    return folders


def _warn_fallback(path: str, what: str) -> None:
    import warnings
    warnings.warn(f"DIFFUSE_PYTHON_PATH={path} has no {what}; using the one installed "
                  "with the diffuse package. Unset DIFFUSE_PYTHON_PATH to silence this.",
                  stacklevel=4)


def find_library() -> Path | None:
    """The C library to use, or None."""
    for folder, packaged in _folders():
        if (folder / LIBRARY).exists():
            if packaged and os.environ.get("DIFFUSE_PYTHON_PATH"):
                _warn_fallback(os.environ["DIFFUSE_PYTHON_PATH"], LIBRARY)
            return folder / LIBRARY
    return None


def load() -> "CtypesCore":
    """The DISCUS core (one per process: the library holds DISCUS's state)."""
    lib = find_library()
    if lib is None:
        raise ImportError(_not_found(LIBRARY))
    return CtypesCore(lib)


def _not_found(what: str) -> str:
    here = Path(__file__).resolve().parent
    path = os.environ.get("DIFFUSE_PYTHON_PATH")
    lines = [
        f"The compiled DISCUS core ({what}) was not found.",
        f"  diffuse package used: {here}",
        f"  DIFFUSE_PYTHON_PATH:  {path or '(not set)'}",
        f"  PYTHONPATH:           {os.environ.get('PYTHONPATH') or '(not set)'}",
        "Install the package with 'pip install .' (in the DiffuseCode folder, in "
        "your virtual environment), or build with -DDIFFUSE_PYTHON=ON and set "
        "DIFFUSE_PYTHON_PATH=<build>/python.",
    ]
    if "site-packages" not in str(here):
        lines.append("The diffuse package above is the source folder, not an installed "
                     "copy: if you installed with pip, unset PYTHONPATH and "
                     "DIFFUSE_PYTHON_PATH (check ~/.bashrc).")
    return "\n".join(lines)


# --- the C library through ctypes ------------------------------------------------
_i32, _i64, _f64, _chars = ctypes.c_int, ctypes.c_int64, ctypes.c_double, ctypes.c_char_p
_P = ctypes.POINTER
_f32a = np.ctypeslib.ndpointer(np.float32, flags="C_CONTIGUOUS")
_i32a = np.ctypeslib.ndpointer(np.int32, flags="C_CONTIGUOUS")

# name: (argument types, result type), as in diffuse_core.h
_SIGNATURES = {
    "dc_api_version": ([], _i32),
    "dc_init": ([], None),
    "dc_reset": ([], None),
    "dc_flush_output": ([], None),
    "dc_run": ([_chars, _i64], _i64),
    "dc_value": ([_chars, _i64], _f64),
    "dc_set_value": ([_chars, _i64, _f64], None),
    "dc_clear_error": ([], None),
    "dc_last_error": ([_P(_i32), _P(_i32), _chars, _chars, _chars], None),
    "dc_seed": ([_i32], None),
    "dc_random_size": ([], _i64),
    "dc_get_random_state": ([_i32a, _i64], _i32),
    "dc_set_random_state": ([_i32a, _i64], _i32),
    "dc_data_sets": ([], _i32),
    "dc_data_length": ([_i32], _i64),
    "dc_get_data": ([_i32, _f32a, _f32a, _i64], _i32),
    "dc_get_data_2d": ([_i32, _f32a, _f32a, _f32a, _i64, _P(_i64), _P(_i64)], _i32),
}


class CtypesCore:
    """DISCUS through the C interface (libdiffuse_core). ctypes releases the
    GIL during every call, so a thread can show the output meanwhile."""

    name = "ctypes"

    def __init__(self, path: Path):
        self.path = Path(path)
        lib = ctypes.CDLL(str(self.path))
        for fname, (args, result) in _SIGNATURES.items():
            f = getattr(lib, fname)
            f.argtypes = args
            f.restype = result
        version = lib.dc_api_version()
        if version != API_VERSION:
            raise ImportError(f"{self.path} has C interface version {version}; this diffuse "
                              f"package needs {API_VERSION}. Rebuild or reinstall.")
        self._lib = lib

    def init(self) -> None:
        self._lib.dc_init()

    def reset(self) -> None:
        self._lib.dc_reset()

    def flush(self) -> None:
        self._lib.dc_flush_output()

    def run(self, text: str) -> None:
        b = text.encode()
        longest = self._lib.dc_run(b, len(b))
        if longest:
            raise ValueError(f"a command line has {longest} characters, more than "
                             "DISCUS's limit of 1024")

    def value(self, expression: str) -> float:
        b = expression.encode()
        return float(self._lib.dc_value(b, len(b)))

    def clear_error(self) -> None:
        self._lib.dc_clear_error()

    def last_error(self) -> tuple[int, str, str, str]:
        """(number, type name, standard message, extra lines as one string
        of 80-character pieces)"""
        num, typ = _i32(), _i32()
        cls, text, msgs = (ctypes.create_string_buffer(n) for n in ERR_LENGTHS)
        self._lib.dc_last_error(ctypes.byref(num), ctypes.byref(typ), cls, text, msgs)
        return (num.value, cls.raw.decode(errors="replace"),
                text.raw.decode(errors="replace"), msgs.raw.decode(errors="replace"))

    def seed(self, value: int) -> None:
        self._lib.dc_seed(int(value))

    def random_size(self) -> int:
        return int(self._lib.dc_random_size())

    def random_state(self) -> np.ndarray:
        n = self.random_size()
        state = np.zeros(n, np.int32)
        if self._lib.dc_get_random_state(state, n):
            raise RuntimeError("dc_get_random_state failed")
        return state

    def set_random_state(self, state: np.ndarray) -> None:
        state = np.ascontiguousarray(state, dtype=np.int32)
        if self._lib.dc_set_random_state(state, state.size):
            raise ValueError(f"state must have {self.random_size()} integers")

    def data_sets(self) -> int:
        return int(self._lib.dc_data_sets())

    def data_length(self, ik: int) -> int:
        return int(self._lib.dc_data_length(int(ik)))

    def data(self, ik: int) -> tuple[np.ndarray, np.ndarray]:
        n = self.data_length(ik)
        x, y = np.zeros(max(n, 0), np.float32), np.zeros(max(n, 0), np.float32)
        if n < 0 or self._lib.dc_get_data(int(ik), x, y, n):
            raise IndexError(f"KUPLOT data set {ik} does not exist")
        return x, y

    def data_2d(self, ik: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, int, int]:
        """x, y, z per point (y running fastest) and the grid size nx, ny"""
        n = self.data_length(ik)
        x, y, z = (np.zeros(max(n, 0), np.float32) for _ in range(3))
        nx, ny = _i64(), _i64()
        if n < 0 or self._lib.dc_get_data_2d(int(ik), x, y, z, n,
                                             ctypes.byref(nx), ctypes.byref(ny)):
            raise ValueError(f"KUPLOT data set {ik} is not a grid")
        return x, y, z, nx.value, ny.value
