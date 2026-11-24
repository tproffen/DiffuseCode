"""Pair distribution function G(r)."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ._session import session
from .crystal import Crystal, _ensure_loaded

_RADIATION = {"xray": "xray", "x-ray": "xray", "neutron": "neutron", "electron": "electron"}


@dataclass
class PDFPattern:
    """A calculated pair distribution function. See :func:`pdf`.

    Attributes
    ----------
    r : numpy.ndarray
        Distance r in Å.
    G : numpy.ndarray
        G(r) = 4πr(ρ(r) − ρ₀) in Å⁻².
    settings : dict
        The arguments it was calculated with.
    crystal : Crystal or None
        The structure it was calculated for.
    """

    r: np.ndarray
    G: np.ndarray
    settings: dict = field(default_factory=dict)
    crystal: Crystal | None = None

    def plot(self, ax=None, **kwargs):
        """Plot G(r) with matplotlib.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes to draw into. Default: a new figure.
        **kwargs
            Passed on to :meth:`matplotlib.axes.Axes.plot`.

        Returns
        -------
        matplotlib.axes.Axes
        """
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots(figsize=(8, 4.5))
        kwargs.setdefault("lw", 1.0)
        kwargs.setdefault("label", (Path(self.crystal.label).stem if self.crystal else None))
        ax.plot(self.r, self.G, **kwargs)
        ax.axhline(0, color="0.6", lw=0.5)
        ax.set_xlabel(r"$r$ (Å)")
        ax.set_ylabel(r"$G(r)$ (Å$^{-2}$)")
        ax.set_xlim(self.r[0], self.r[-1])
        return ax


def pdf(crystal: Crystal, *, rmax: float = 20.0, dr: float = 0.01,
        radiation: str = "xray", qmax: float = 0.0, qdamp: float = 0.0,
        qbroad: float = 0.0, thermal: bool = True, corrlin: float = 0.0,
        corrquad: float = 0.0, particle_diameter: float | None = None) -> PDFPattern:
    """Calculate the pair distribution function G(r) of a crystal with DISCUS.

    Parameters
    ----------
    crystal : Crystal
        The structure, from :func:`read_structure`.
    rmax, dr : float
        Largest distance and step in Å.
    radiation : {"xray", "neutron", "electron"}
        Weights the atom pairs by scattering power.
    qmax : float
        Qmax in Å⁻¹ of the termination ripples to include; 0 means none.
    qdamp, qbroad : float
        Instrumental damping and resolution broadening (as in PDFgui).
    thermal : bool
        Broaden each pair by a Gaussian from the displacement parameters.
        False gives sharp, histogram-like peaks.
    corrlin, corrquad : float
        Linear and quadratic correlated-motion corrections (peak sharpening
        at small r); only with ``thermal=True``.
    particle_diameter : float, optional
        Treat the crystal as a finite spherical particle of this diameter in
        Å (no periodic boundaries; the -4πρ₀r baseline is corrected for a
        sphere). Default: an infinite crystal with periodic boundaries.

    Returns
    -------
    PDFPattern

    Raises
    ------
    ValueError
        For an unknown ``radiation`` or invalid ranges.
    MacroError
        If the calculation fails inside DISCUS.

    Examples
    --------
    >>> g = dc.pdf(cr, rmax=20, dr=0.01, qmax=25)  # doctest: +SKIP
    >>> g.plot()  # doctest: +SKIP
    """
    try:
        rad = _RADIATION[radiation.lower()]
    except KeyError:
        raise ValueError(f"unknown radiation '{radiation}'; use {sorted(set(_RADIATION))}") from None
    if not 0 < dr < rmax:
        raise ValueError("need 0 < dr < rmax")
    if qmax < 0 or qdamp < 0 or qbroad < 0:
        raise ValueError("qmax, qdamp and qbroad must not be negative")
    if particle_diameter is not None and particle_diameter <= 0:
        raise ValueError("particle_diameter must be positive")

    if particle_diameter is None:
        boundary = "set bound, periodic, 3D\n  set finite, periodic"
    else:
        boundary = f"set bound, crystal\n  set finite, sphere, {particle_diameter}"

    _ensure_loaded(crystal)
    s = session()
    out = s.workdir / "pdf.gr"
    out.unlink(missing_ok=True)
    s.run(f"""
discus
pdf
  reset all
  set range, {rmax}, {dr}
  set qmax, {qmax}
  set qdamp, {qdamp}
  set qbroad, {qbroad}
  {boundary}
  set density, auto
  set corrlin, {corrlin}
  set corrquad, {corrquad}
  set therm, {"gauss" if thermal else "off"}
  set weight, 1.0
  set srat, 1.0, 2.5
  set rad, {rad}
  calc
  save pdf, {out.name}
exit
exit
""", cwd=s.workdir)
    data = np.loadtxt(out, comments="#", ndmin=2)
    settings = dict(rmax=rmax, dr=dr, radiation=rad, qmax=qmax, qdamp=qdamp,
                    qbroad=qbroad, thermal=thermal, corrlin=corrlin,
                    corrquad=corrquad, particle_diameter=particle_diameter)
    return PDFPattern(r=data[:, 0], G=data[:, 1], settings=settings, crystal=crystal)
