"""
Pair distribution functions
===========================

Calculate the pair distribution function G(r) with :func:`diffuse.pdf`:
how a finite Q range adds termination ripples, and how X-rays and neutrons
weight the atom pairs differently.

To run it yourself, from the repository root::

    DIFFUSE_PYTHON_PATH=build/python PYTHONPATH=python python3 python/examples/pdf_demo.py
"""

from pathlib import Path

import matplotlib.pyplot as plt

import diffuse as dc

# sphinx_gallery_start_ignore
try:
    HERE = Path(__file__).resolve().parent
except NameError:            # the docs gallery runs examples from their folder
    HERE = Path.cwd()
# sphinx_gallery_end_ignore

# %%
# Nickel: the effect of a finite Q range
# --------------------------------------
# ``qmax=0`` gives the PDF of an ideal measurement. A real one stops at
# some Qmax, which adds ripples around each peak; ``qdamp`` describes the
# instrument resolution, which damps the peaks at large r.

nickel = dc.read_structure(HERE / "nickel.cell")
ideal = dc.pdf(nickel, rmax=20, dr=0.01)
measured = dc.pdf(nickel, rmax=20, dr=0.01, qmax=20, qdamp=0.04)

ax = ideal.plot(label="ideal (no Qmax)")
measured.plot(ax=ax, label="Qmax = 20 Å⁻¹, Qdamp = 0.04 Å⁻¹")
ax.set_title("Nickel (fcc)")
ax.legend()

# %%
# The first peak sits at the nearest-neighbour distance a/√2:

first = ideal.r[(ideal.G * (ideal.r < 3)).argmax()]
print(f"first peak at r = {first:.3f} Å, a/√2 = {nickel.a / 2**0.5:.3f} Å")

# %%
# Anatase: X-rays and neutrons
# ----------------------------
# Each atom pair is weighted by the scattering power of its two atoms.
# Titanium dominates for X-rays; for neutrons its scattering length is
# negative, so Ti–O peaks point down.

anatase = dc.read_structure(HERE / "anatase.cif")
xray = dc.pdf(anatase, rmax=10, dr=0.01, qmax=25, radiation="xray")
neutron = dc.pdf(anatase, rmax=10, dr=0.01, qmax=25, radiation="neutron")

fig, ax = plt.subplots(figsize=(8, 4.5))
xray.plot(ax=ax, label="X-ray")
neutron.plot(ax=ax, label="neutron")
ax.set_title("Anatase (TiO₂)")
ax.legend()

plt.show()
