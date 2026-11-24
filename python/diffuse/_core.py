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
API_VERSION = 4              # DC_API_VERSION in diffuse_core.h
ERR_LENGTHS = (4, 45, 560)   # DC_ERR_CLASS_LEN, DC_ERR_TEXT_LEN, DC_ERR_MSGS_LEN
SPCGR_LEN = 16               # DC_SPCGR_LEN
NAME_LEN = 4                 # DC_NAME_LEN


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
_i64a = np.ctypeslib.ndpointer(np.int64, flags="C_CONTIGUOUS")
_f64a = np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS")

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
    "dc_crystal_size": ([_P(_i64), _P(_i32)], None),
    "dc_get_cell": ([_f64a, _chars, _i64a, _P(_i64)], None),
    "dc_get_atoms": ([_i64, _f64a, _i32a, _i32a, _i32a], _i32),
    "dc_get_atom_types": ([_i32, _chars, _f64a, _f64a], _i32),
    "dc_set_crystal": ([_f64a, _chars, _i64, _i64a, _i64, _i32, _chars, _f64a, _f64a,
                        _i32, _f64a, _i64, _f64a, _i32a, _i32a], _i32),
    "dc_get_type_anis": ([_i32, _f64a], _i32),
    "dc_adp_count": ([], _i64),
    "dc_get_adps": ([_i64, _f64a, _f64a], _i32),
    "dc_get_atom_adps": ([_i64, _i32a], _i32),
    "dc_molecule_size": ([_P(_i64), _P(_i32), _P(_i64)], None),
    "dc_get_molecules": ([_i64, _i64, _i32a, _i32a, _i32a, _i32a], _i32),
    "dc_set_molecules": ([_i64, _i64, _i32a, _i32a, _i32a], _i32),
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

    # -- the crystal DISCUS holds ---------------------------------------------
    def crystal_size(self) -> tuple[int, int]:
        """(number of atoms, number of atom types); types run 0..ntypes."""
        natoms, ntypes = _i64(), _i32()
        self._lib.dc_crystal_size(ctypes.byref(natoms), ctypes.byref(ntypes))
        return natoms.value, ntypes.value

    def cell(self) -> tuple[tuple[float, ...], str, tuple[int, int, int], int]:
        """(a, b, c, alpha, beta, gamma), space group, unit cells, atoms per cell"""
        lattice = np.zeros(6, np.float64)
        spcgr = ctypes.create_string_buffer(SPCGR_LEN)
        icc = np.zeros(3, np.int64)
        ncatoms = _i64()
        self._lib.dc_get_cell(lattice, spcgr, icc, ctypes.byref(ncatoms))
        return (tuple(float(v) for v in lattice), spcgr.raw.decode(errors="replace").strip(),
                tuple(int(v) for v in icc), ncatoms.value)

    def atoms(self) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Positions (natoms, 3), atom types, property flags, molecule numbers"""
        n = self.crystal_size()[0]
        pos = np.zeros((n, 3), np.float64)
        types, prop, mole = (np.zeros(n, np.int32) for _ in range(3))
        if self._lib.dc_get_atoms(n, pos, types, prop, mole):
            raise RuntimeError("dc_get_atoms: the number of atoms changed")
        return pos, types, prop, mole

    def atom_types(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Names, B (A^2) and occupancy of atom types 0..ntypes"""
        ntypes = self.crystal_size()[1]
        names = ctypes.create_string_buffer(NAME_LEN * (ntypes + 1))
        biso, occ = np.zeros(ntypes + 1, np.float64), np.zeros(ntypes + 1, np.float64)
        if self._lib.dc_get_atom_types(ntypes, names, biso, occ):
            raise RuntimeError("dc_get_atom_types: the number of atom types changed")
        raw = names.raw.decode(errors="replace")
        labels = np.array([raw[i:i + NAME_LEN].strip() for i in range(0, len(raw), NAME_LEN)])
        return labels, biso, occ

    def set_crystal(self, lattice, space_group: str, cells, atoms_per_cell: int,
                    names, b_iso, occupancy, positions, types, properties,
                    u_aniso=None) -> bool:
        """Replace DISCUS's crystal. names/b_iso/occupancy (and u_aniso,
        shape (ntypes, 6), or None) per type 1..ntypes, positions (natoms, 3),
        types 1..ntypes. False on a DISCUS error."""
        names = list(names)
        if any(len(n.encode()) > NAME_LEN for n in names):
            raise ValueError(f"atom names have at most {NAME_LEN} characters: {names}")
        block = b"".join(n.encode().ljust(NAME_LEN) for n in names)
        sg = space_group.encode()
        pos = np.ascontiguousarray(positions, dtype=np.float64).reshape(-1, 3)
        return self._lib.dc_set_crystal(
            np.ascontiguousarray(lattice, dtype=np.float64), sg, len(sg),
            np.ascontiguousarray(cells, dtype=np.int64), int(atoms_per_cell),
            len(names), block, np.ascontiguousarray(b_iso, dtype=np.float64),
            np.ascontiguousarray(occupancy, dtype=np.float64),
            0 if u_aniso is None else 1,
            np.ascontiguousarray(np.zeros((len(names), 6)) if u_aniso is None else u_aniso,
                                 dtype=np.float64).reshape(-1),
            len(pos), pos,
            np.ascontiguousarray(types, dtype=np.int32),
            np.ascontiguousarray(properties, dtype=np.int32)) == 0

    # -- displacement parameters ---------------------------------------------
    def type_anis(self) -> np.ndarray:
        """ADPs of the types 0..ntypes as given, shape (ntypes+1, 6)"""
        ntypes = self.crystal_size()[1]
        out = np.zeros((ntypes + 1, 6), np.float64)
        if self._lib.dc_get_type_anis(ntypes, out):
            raise RuntimeError("dc_get_type_anis: the number of atom types changed")
        return out

    def adps(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Distinct tensors (nanis, 6), principal axes (nanis, 3, 4), and per
        atom the tensor (1-based index)"""
        nanis = int(self._lib.dc_adp_count())
        uij, prin = np.zeros((nanis, 6)), np.zeros((nanis, 3, 4))
        if self._lib.dc_get_adps(nanis, uij, prin):
            raise RuntimeError("dc_get_adps: the number of tensors changed")
        n = self.crystal_size()[0]
        index = np.zeros(n, np.int32)
        if self._lib.dc_get_atom_adps(n, index):
            raise RuntimeError("dc_get_atom_adps: the number of atoms changed")
        return uij, prin, index

    # -- molecules ---------------------------------------------------------------
    def molecules(self) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Per molecule: type, character, number of atoms; and all their atoms
        (1-based), one molecule after the other"""
        nmole, ntypes, natoms = _i64(), _i32(), _i64()
        self._lib.dc_molecule_size(ctypes.byref(nmole), ctypes.byref(ntypes), ctypes.byref(natoms))
        mtype, mchar, mlen = (np.zeros(nmole.value, np.int32) for _ in range(3))
        atoms = np.zeros(natoms.value, np.int32)
        if self._lib.dc_get_molecules(nmole.value, natoms.value, mtype, mchar, mlen, atoms):
            raise RuntimeError("dc_get_molecules: sizes changed")
        return mtype, mchar, mlen, atoms

    def set_molecules(self, types, lengths, atoms) -> bool:
        """Replace the molecules: per molecule its type (1, 2 ...) and number
        of atoms, and all their atoms (1-based). False on a DISCUS error."""
        types = np.ascontiguousarray(types, dtype=np.int32)
        lengths = np.ascontiguousarray(lengths, dtype=np.int32)
        atoms = np.ascontiguousarray(atoms, dtype=np.int32)
        return self._lib.dc_set_molecules(len(types), len(atoms), types, lengths, atoms) == 0
