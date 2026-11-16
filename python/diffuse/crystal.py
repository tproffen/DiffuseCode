"""Reading structures into DISCUS."""

from __future__ import annotations

import os
import shutil
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ._session import DiffuseError, session

# Files that DISCUS reads as an asymmetric unit (then applies the space
# group and builds the supercell) and files it reads as a whole crystal.
_CELL_SUFFIXES = {".cell", ".cif"}
_STRU_SUFFIXES = {".stru"}


@dataclass(eq=False)
class Crystal:
    """A structure loaded into DISCUS.

    DISCUS holds one structure at a time. A Crystal remembers how it was
    made, so a calculation on an older Crystal reads it in again first.
    """

    path: Path
    supercell: tuple[int, int, int] | None
    natoms: int
    lattice: tuple[float, float, float, float, float, float]  # a, b, c, alpha, beta, gamma
    _atoms: tuple[np.ndarray, np.ndarray] | None = field(default=None, repr=False)

    @property
    def a(self) -> float:
        return self.lattice[0]

    @property
    def volume(self) -> float:
        """Volume of one unit cell in Å^3."""
        a, b, c, al, be, ga = self.lattice
        ca, cb, cg = np.cos(np.radians([al, be, ga]))
        return float(a * b * c * np.sqrt(1 - ca**2 - cb**2 - cg**2 + 2 * ca * cb * cg))

    def atoms(self) -> tuple[np.ndarray, np.ndarray]:
        """Atom names and fractional coordinates, shape (natoms,) and (natoms, 3).

        DISCUS saves the crystal to a ``.stru`` file, which is read back here.
        """
        if self._atoms is None:
            _ensure_loaded(self)
            s = session()
            out = s.workdir / "atoms.stru"
            out.unlink(missing_ok=True)
            s.run(f"discus\n  save\n    outfile {out.name}\n    write all\n    run\n  exit\nexit",
                  cwd=s.workdir)
            lines = out.read_text().splitlines()
            start = next(i for i, line in enumerate(lines) if line.split()[:1] == ["atoms"]) + 1
            rows = [line.replace(",", " ").split() for line in lines[start:] if line.strip()]
            names = np.array([r[0] for r in rows])
            xyz = np.array([[float(v) for v in r[1:4]] for r in rows])
            if len(names) != self.natoms:
                raise DiffuseError(f"saved structure has {len(names)} atoms, expected {self.natoms}")
            self._atoms = (names, xyz)
        return self._atoms

    def positions(self) -> np.ndarray:
        """Fractional coordinates, shape (natoms, 3)."""
        return self.atoms()[1]

    def __repr__(self) -> str:
        cell = ", ".join(f"{v:g}" for v in self.lattice)
        sc = f", supercell={self.supercell}" if self.supercell else ""
        return f"Crystal('{self.path.name}'{sc}, natoms={self.natoms}, lattice=({cell}))"


def read_structure(path: str | os.PathLike,
                   supercell: tuple[int, int, int] | int | None = None) -> Crystal:
    """Read a structure file into DISCUS.

    ``.cell`` and ``.cif`` files hold an asymmetric unit: DISCUS applies the
    space group and builds ``supercell`` unit cells (default 1 x 1 x 1).
    ``.stru`` files hold a whole crystal and are read as they are, unless a
    ``supercell`` is given.
    """
    path = Path(path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(path)
    suffix = path.suffix.lower()
    if suffix not in _CELL_SUFFIXES | _STRU_SUFFIXES:
        raise DiffuseError(f"unsupported structure file '{path.name}' "
                           f"(use {', '.join(sorted(_CELL_SUFFIXES | _STRU_SUFFIXES))})")
    if isinstance(supercell, int):
        supercell = (supercell,) * 3
    if supercell is not None:
        supercell = tuple(int(n) for n in supercell)
        if len(supercell) != 3 or min(supercell) < 1:
            raise ValueError(f"supercell must be three positive integers, not {supercell}")

    if suffix in _STRU_SUFFIXES and supercell is None:
        command = f"stru {path.name}"
    else:
        nx, ny, nz = supercell or (1, 1, 1)
        command = f"cell {path.name}, {nx}, {ny}, {nz}"

    # DISCUS may write files next to the input (a CIF is converted to .stru
    # first), so read a copy in the session's scratch folder.
    s = session()
    inputs = s.workdir / "input"
    shutil.rmtree(inputs, ignore_errors=True)
    inputs.mkdir()
    shutil.copy(path, inputs / ("structure" + suffix))
    command = command.replace(path.name, "structure" + suffix, 1)

    # n[] and lat[] exist only inside the discus section: copy them to the
    # r[] entries reserved for this package (r[490] ... r[498]).
    s.run(f"""
discus
  read
    {command}
  r[490] = n[1]
  r[491] = lat[1]
  r[492] = lat[2]
  r[493] = lat[3]
  r[494] = lat[4]
  r[495] = lat[5]
  r[496] = lat[6]
exit
""", cwd=inputs)
    crystal = Crystal(
        path=path,
        supercell=supercell,
        natoms=int(s.value("r[490]")),
        lattice=tuple(s.value(f"r[{i}]") for i in range(491, 497)),
    )
    s.current_crystal = crystal
    return crystal


def _ensure_loaded(crystal: Crystal) -> None:
    """Make ``crystal`` the structure DISCUS holds, reading it again if needed."""
    s = session()
    if s.current_crystal is crystal:
        return
    fresh = read_structure(crystal.path, crystal.supercell)
    if fresh.natoms != crystal.natoms:
        raise DiffuseError(f"{crystal.path.name} now gives {fresh.natoms} atoms, "
                           f"not {crystal.natoms}; has the file changed?")
    s.current_crystal = crystal
