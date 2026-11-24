"""
Crystals from Python: arrays and ASE
====================================

DISCUS can take a crystal straight from Python arrays: positions, names,
displacement parameters. That makes it easy to build or change a structure
with NumPy (or with ASE, the Atomic Simulation Environment) and calculate
its scattering with DISCUS.

To run it yourself, from the repository root::

    DIFFUSE_PYTHON_PATH=build/python PYTHONPATH=python python3 python/examples/crystals_from_python.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from ase.cluster import Octahedron

import diffuse as dc

# sphinx_gallery_start_ignore
try:
    HERE = Path(__file__).resolve().parent
except NameError:            # the docs gallery runs examples from their folder
    HERE = Path.cwd()
# sphinx_gallery_end_ignore

# %%
# A crystal from NumPy arrays
# ---------------------------
# Copper, face-centred cubic, 6 x 6 x 6 unit cells. ``positions`` are
# fractional coordinates in units of the unit cell, so they run from 0 to 6.

fcc = np.array([[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
cells = np.array([[i, j, k] for i in range(6) for j in range(6) for k in range(6)])
xyz = (cells[:, None, :] + fcc[None, :, :]).reshape(-1, 3)
copper = dc.make_crystal((3.615, 3.615, 3.615, 90, 90, 90), ["CU"] * len(xyz), xyz,
                         space_group="F m -3 m", b_iso=0.55, cells=(6, 6, 6))
print(copper)

# %%
# Changing a structure in Python
# ------------------------------
# Read anatase, take its atoms out as arrays, remove 10 % of the oxygen at
# random and give the result back to DISCUS. The powder patterns show what
# the vacancies do.

anatase = dc.read_structure(HERE / "anatase.cif", supercell=(8, 8, 3))
names, xyz = anatase.atoms()
rng = np.random.default_rng(1)
oxygen = np.flatnonzero(names == "O2-")
keep = np.ones(len(names), bool)
keep[rng.choice(oxygen, size=len(oxygen) // 10, replace=False)] = False
b = dict(zip(*anatase.atom_types()[:2]))
defective = dc.make_crystal(anatase.lattice, names[keep], xyz[keep], b_iso=b,
                            space_group=anatase.space_group)
print(anatase.natoms, "atoms ->", defective.natoms, "atoms")

kw = dict(qmin=0.5, qmax=8, dq=0.01, radiation="neutron")
perfect, vacancies = dc.powder(anatase, **kw), dc.powder(defective, **kw)

fig, ax = plt.subplots(figsize=(7, 3.5))
ax.plot(perfect.q, perfect.intensity, label="anatase")
ax.plot(vacancies.q, vacancies.intensity * perfect.intensity.max() / vacancies.intensity.max(),
        label="10 % O vacancies (scaled)", lw=1)
ax.set_xlabel("Q (Å$^{-1}$)")
ax.set_ylabel("intensity (neutrons)")
ax.legend()
fig.tight_layout()

# %%
# A nanoparticle from ASE
# -----------------------
# ASE builds the particle (a gold octahedron); :func:`diffuse.from_ase`
# hands it to DISCUS, which calculates its PDF. ``Crystal.to_ase()`` goes
# the other way.

particle = Octahedron("Au", length=7, cutoff=2, latticeconstant=4.078)
particle.center(vacuum=5.0)                 # a box around the particle
gold = dc.from_ase(particle, b_iso=0.6)
print(gold.natoms, "gold atoms in a", gold.space_group, "box")

g = dc.pdf(gold, rmax=20, dr=0.02, radiation="xray")
fig, ax = plt.subplots(figsize=(7, 3.5))
ax.plot(g.r, g.G)
ax.set_xlabel("r (Å)")
ax.set_ylabel("G(r)")
ax.set_title("Au octahedron")
fig.tight_layout()

back = gold.to_ase()
print("back in ASE:", back.get_chemical_formula())
