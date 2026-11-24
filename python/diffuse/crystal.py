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
    """A structure loaded into DISCUS. Create it with :func:`read_structure`.

    DISCUS holds one structure at a time. A ``Crystal`` remembers how it was
    made, so a calculation on an older ``Crystal`` reads it in again first.

    Attributes
    ----------
    path : pathlib.Path
        The structure file it was read from.
    supercell : tuple of int or None
        Number of unit cells along a, b and c, or None for a ``.stru`` file
        read as it is.
    natoms : int
        Number of atoms in the whole crystal.
    lattice : tuple of float
        ``(a, b, c, alpha, beta, gamma)`` in Å and degrees.
    space_group : str
        Space group symbol as DISCUS uses it, for example ``"F m -3 m"``.
    """

    path: Path
    supercell: tuple[int, int, int] | None
    natoms: int
    lattice: tuple[float, float, float, float, float, float]  # a, b, c, alpha, beta, gamma
    space_group: str = "P1"
    _atoms: tuple[np.ndarray, np.ndarray] | None = field(default=None, repr=False)
    _types: tuple[np.ndarray, np.ndarray, np.ndarray] | None = field(default=None, repr=False)

    @property
    def a(self) -> float:
        """Lattice parameter a in Å."""
        return self.lattice[0]

    @property
    def volume(self) -> float:
        """Volume of one unit cell in Å^3."""
        a, b, c, al, be, ga = self.lattice
        ca, cb, cg = np.cos(np.radians([al, be, ga]))
        return float(a * b * c * np.sqrt(1 - ca**2 - cb**2 - cg**2 + 2 * ca * cb * cg))

    def atoms(self) -> tuple[np.ndarray, np.ndarray]:
        """Atom names and fractional coordinates.

        Returns
        -------
        names : numpy.ndarray of str, shape (natoms,)
            Atom names as DISCUS uses them, for example ``"TI4+"``.
        xyz : numpy.ndarray of float, shape (natoms, 3)
            Fractional coordinates in units of the unit cell.

        Notes
        -----
        Read straight from DISCUS's arrays (no file). The result is cached
        on the ``Crystal``.
        """
        if self._atoms is None:
            self._read_arrays()
        return self._atoms

    def atom_types(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """The atom types: names, isotropic displacement parameters and occupancies.

        Returns
        -------
        names : numpy.ndarray of str, shape (ntypes,)
            Atom type names, for example ``"TI4+"``, in DISCUS's order (the
            names :meth:`atoms` uses).
        b_iso : numpy.ndarray of float, shape (ntypes,)
            Isotropic displacement parameter B in Å^2 (B = 8 pi^2 <u^2>).
        occupancy : numpy.ndarray of float, shape (ntypes,)
            Occupancy of each type.

        Notes
        -----
        Vacancies (DISCUS's type 0, ``"VOID"``) are left out.
        """
        if self._types is None:
            self._read_arrays()
        return self._types

    def _read_arrays(self) -> None:
        """Fill the caches from DISCUS's crystal arrays (C interface)."""
        _ensure_loaded(self)
        core = session()._c
        pos, types, _prop, _mole = core.atoms()
        names, b_iso, occupancy = core.atom_types()
        if len(pos) != self.natoms:
            raise DiffuseError(f"DISCUS holds {len(pos)} atoms, expected {self.natoms}")
        self._atoms = (names[types], pos)
        self._types = (names[1:], b_iso[1:], occupancy[1:])

    def positions(self) -> np.ndarray:
        """Fractional coordinates, shape (natoms, 3). See :meth:`atoms`."""
        return self.atoms()[1]

    def lattice_vectors(self) -> np.ndarray:
        """The unit-cell vectors a, b, c as rows, in Å, shape (3, 3).

        Standard orientation: a along x, b in the xy plane, c completing a
        right-handed set.
        """
        a, b, c, al, be, ga = self.lattice
        ca, cb, cg = np.cos(np.radians([al, be, ga]))
        sg = np.sin(np.radians(ga))
        cx = c * cb
        cy = c * (ca - cb * cg) / sg
        cz = np.sqrt(max(c * c - cx * cx - cy * cy, 0.0))
        return np.array([[a, 0.0, 0.0], [b * cg, b * sg, 0.0], [cx, cy, cz]])

    def cartesian(self) -> np.ndarray:
        """Atom positions in Å, shape (natoms, 3), in the orientation of
        :meth:`lattice_vectors`."""
        return self.positions() @ self.lattice_vectors()

    def view(self, style: str = "ball", *, scale: float = 1.0, cell: str | None = "crystal",
             show=None, colors: dict | None = None, width: int = 640, height: int = 480,
             background: str = "white"):
        """Show the structure in 3D in a Jupyter notebook (py3Dmol).

        Rotate with the left mouse button, zoom with the wheel, move with the
        right button (or Ctrl + left). Uses ``py3Dmol`` (3Dmol.js), which
        is installed with the package.

        Parameters
        ----------
        style : {"ball", "spacefill", "stick", "ballstick"}
            Small spheres (default), touching spheres, or bonds; bonds are
            drawn between atoms closer than the sum of their covalent radii.
        scale : float
            Size factor for the spheres and sticks, e.g. 3 for a sparse
            model structure whose atoms are far apart.
        cell : {"crystal", "unit"} or None
            Draw the edges of the whole crystal (all unit cells), of one
            unit cell, or nothing.
        show : str or sequence of str, optional
            Only these atom types (DISCUS names such as ``"CU"``), for
            example to see how one kind of atom is ordered.
        colors : dict, optional
            Colours for atom types, e.g. ``{"CU": "orange", "AU": "gold"}``.
            Other atoms are coloured by element (Jmol colours).
        width, height : int
            Size in pixels.
        background : str
            Background colour.

        Returns
        -------
        StructureView
            Displayed when it is the last expression in a notebook cell
            (or call ``.show()``). ``.write_html("file.html")`` saves a page
            that shows the structure in any browser. The py3Dmol view is
            ``.viewer``; other py3Dmol methods (``addLabel``, ``setView`` ...)
            can be called on the StructureView directly, then ``.show()``.

        Notes
        -----
        Fine for tens of thousands of atoms with ``"ball"``; bonds
        (``"stick"``, ``"ballstick"``) get slow for large crystals.
        """
        try:
            import py3Dmol
        except ImportError:
            raise ImportError("Crystal.view() needs py3Dmol, which is installed with "
                              "the diffuse package: pip install py3Dmol") from None
        if style not in _VIEW_STYLES:
            raise ValueError(f"style must be one of {sorted(_VIEW_STYLES)}, not {style!r}")
        if not scale > 0:
            raise ValueError("scale must be positive")
        look = {kind: {k: v * scale for k, v in spec.items()}
                for kind, spec in _VIEW_STYLES[style].items()}
        if cell not in ("crystal", "unit", None):
            raise ValueError(f"cell must be 'crystal', 'unit' or None, not {cell!r}")
        names, frac = self.atoms()
        keep = np.ones(len(names), bool)
        if show is not None:
            wanted = {show.upper()} if isinstance(show, str) else {s.upper() for s in show}
            unknown = wanted - set(names)
            if unknown:
                raise ValueError(f"no atoms of type {sorted(unknown)}; the crystal has {sorted(set(names))}")
            keep = np.isin(names, list(wanted))
        names = names[keep]
        xyz = frac[keep] @ self.lattice_vectors()

        model = [str(len(names)), self.path.name]
        model += [f"{_element(n)} {x:.5f} {y:.5f} {z:.5f}" for n, (x, y, z) in zip(names, xyz)]
        v = py3Dmol.view(width=width, height=height)
        v.addModel("\n".join(model) + "\n", "xyz")
        v.setStyle({}, look)
        for name, color in (colors or {}).items():
            index = np.nonzero(names == name.upper())[0]
            if len(index):
                spec = {k: dict(s, color=color) for k, s in look.items()}
                v.setStyle({"index": index.tolist()}, spec)
        if cell is not None:
            if cell == "unit":
                lo, n = np.zeros(3), np.ones(3)
            else:          # all unit cells: from floor(min) to floor(max) + 1
                lo = np.floor(frac.min(axis=0) + 1e-6)
                n = np.floor(frac.max(axis=0) + 1e-6) + 1 - lo
            for start, end in _box_edges(lo, n):
                p, q = start @ self.lattice_vectors(), end @ self.lattice_vectors()
                v.addLine({"start": dict(zip("xyz", map(float, p))),
                           "end": dict(zip("xyz", map(float, q))),
                           "color": "black", "dashed": False})
        v.setBackgroundColor(background)
        v.zoomTo()
        return StructureView(v)


class StructureView:
    """A 3D view of a crystal, from :meth:`Crystal.view`.

    Wraps a ``py3Dmol.view`` (in ``viewer``) so that a notebook shows the
    view and nothing else; any other py3Dmol method is passed on to it.
    """

    def __init__(self, viewer):
        self.viewer = viewer

    def show(self):
        """Display the view (in a notebook)."""
        self.viewer.show()

    def write_html(self, path=None, fullpage: bool = True):
        """Save the view as a web page (or return the HTML if ``path`` is
        None). The page loads 3Dmol.js from the internet when opened."""
        return self.viewer.write_html(path, fullpage=fullpage)

    def _ipython_display_(self):
        self.viewer.show()

    def __getattr__(self, name):
        return getattr(self.viewer, name)

    def __repr__(self) -> str:
        cell = ", ".join(f"{v:g}" for v in self.lattice)
        sc = f", supercell={self.supercell}" if self.supercell else ""
        return f"Crystal('{self.path.name}'{sc}, natoms={self.natoms}, lattice=({cell}))"


_VIEW_STYLES = {
    "ball": {"sphere": {"scale": 0.3}},
    "spacefill": {"sphere": {"scale": 1.0}},
    "stick": {"stick": {"radius": 0.15}},
    "ballstick": {"sphere": {"scale": 0.25}, "stick": {"radius": 0.12}},
}


def _element(name: str) -> str:
    """Element symbol from a DISCUS atom name: 'TI4+' -> 'Ti', 'O2-' -> 'O'."""
    letters = ""
    for ch in name:
        if not ch.isalpha():
            break
        letters += ch
    letters = letters[:2] or "X"
    return letters[0].upper() + letters[1:].lower()


def _box_edges(lo, n):
    """The 12 edges of the box lo .. lo + n, in fractional coordinates."""
    corners = [lo + n * np.array(c) for c in np.ndindex(2, 2, 2)]
    return [(corners[i], corners[j]) for i in range(8) for j in range(i + 1, 8)
            if np.count_nonzero(np.array(np.unravel_index(i, (2, 2, 2)))
                                != np.array(np.unravel_index(j, (2, 2, 2)))) == 1]


def read_structure(path: str | os.PathLike,
                   supercell: tuple[int, int, int] | int | None = None) -> Crystal:
    """Read a structure file into DISCUS.

    Parameters
    ----------
    path : str or os.PathLike
        A ``.cif`` or ``.cell`` file (an asymmetric unit: DISCUS applies the
        space group) or a ``.stru`` file (a whole crystal).
    supercell : int or tuple of 3 int, optional
        Number of unit cells along a, b and c; an int means the same along
        all three. Default: 1 x 1 x 1 for ``.cif``/``.cell``, and a ``.stru``
        file is read as it is.

    Returns
    -------
    Crystal
        The structure, now held by DISCUS.

    Raises
    ------
    FileNotFoundError
        If ``path`` does not exist.
    DiffuseError
        If the file type is not supported.
    MacroError
        If DISCUS cannot read the file.

    Examples
    --------
    >>> import diffuse as dc
    >>> cr = dc.read_structure("anatase.cif", supercell=(10, 10, 4))  # doctest: +SKIP
    >>> cr.natoms  # doctest: +SKIP
    4800
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

    s.run(f"discus\n  read\n    {command}\nexit\n", cwd=inputs)
    lattice, space_group, _cells, _per_cell = s._c.cell()
    crystal = Crystal(
        path=path,
        supercell=supercell,
        natoms=s._c.crystal_size()[0],
        lattice=lattice,
        space_group=space_group,
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
