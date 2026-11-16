"""
Single-crystal diffraction and diffuse scattering
=================================================

Calculate intensities on planes in reciprocal space with
:func:`diffuse.fourier`, for an ordered crystal (anatase) and for a Cu/Au
layer with short-range order, where removing the average structure leaves
the diffuse scattering.

To run it yourself, from the repository root::

    DIFFUSE_PYTHON_PATH=build/python PYTHONPATH=python python3 python/examples/fourier_demo.py
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
# Bragg peaks of an ordered crystal
# ---------------------------------
# The hk0 plane of anatase (I4\ :sub:`1`/amd), 8 x 8 x 4 unit cells. The
# corners give the plane in reciprocal-lattice units: ``ll`` lower left,
# ``lr`` lower right, ``ul`` upper left.

anatase = dc.read_structure(HERE / "anatase.cif", supercell=(8, 8, 4))
hk0 = dc.fourier(anatase, ll=(-4, -4, 0), lr=(4, -4, 0), ul=(-4, 4, 0),
                 npoints=(321, 321), radiation="xray")
ax = hk0.plot(log=True, cmap="magma")
ax.set_title("Anatase, hk0 plane (X-rays, log scale)")
ax.set_aspect("equal")

# %%
# The result is a NumPy array plus its reciprocal-space coordinates, so
# you can work with it directly, for example to find the strongest
# reflection apart from the origin:

import numpy as np

hkl = hk0.hkl()                                   # shape (321, 321, 3)
away = np.linalg.norm(hkl, axis=-1) > 0.5         # leave out 000
i, j = np.unravel_index(np.where(away, hk0.intensity, 0).argmax(), hk0.intensity.shape)
print("strongest reflection:", hkl[i, j].round(3))

# %%
# Diffuse scattering of a disordered layer
# ----------------------------------------
# A single layer of 50 x 50 cells, occupied by Cu and about 35 % Au. The
# Au atoms cluster: a neighbour of an Au atom is more likely to be Au than
# the average composition says (short-range order). The total scattering
# is dominated by the Bragg peaks of the average structure;
# ``subtract_average=100`` removes them and leaves the diffuse scattering.
# Because like atoms cluster, it gathers around the Bragg positions.

layer = dc.read_structure(HERE / "cuau_layer.stru")
plane = dict(ll=(0, 0, 0), lr=(3, 0, 0), ul=(0, 3, 0), npoints=(241, 241))
total = dc.fourier(layer, **plane)
diffuse = dc.fourier(layer, subtract_average=100, **plane)

fig, (left, right) = plt.subplots(1, 2, figsize=(11, 4.6), constrained_layout=True)
total.plot(ax=left, log=True, cmap="magma")
left.set_title("Total scattering (log scale)")
diffuse.plot(ax=right, cmap="magma")
right.set_title("Diffuse scattering (average removed)")
for ax in (left, right):
    ax.set_aspect("equal")

plt.show()
