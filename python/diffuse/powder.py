"""Powder diffraction patterns."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ._session import session
from .crystal import Crystal, _ensure_loaded

_RADIATION = {"xray": "xray", "x-ray": "xray", "neutron": "neutron", "electron": "electron"}
_METHOD = {"debye": "debye", "complete": "compl"}


@dataclass
class PowderPattern:
    """A calculated powder pattern on an evenly spaced Q grid. See :func:`powder`.

    Attributes
    ----------
    q : numpy.ndarray
        Momentum transfer Q in Å⁻¹.
    intensity : numpy.ndarray
        Intensity at each Q, arbitrary units.
    settings : dict
        The arguments the pattern was calculated with.
    crystal : Crystal or None
        The structure it was calculated for.
    """

    q: np.ndarray                   # Å^-1
    intensity: np.ndarray
    settings: dict = field(default_factory=dict)
    crystal: Crystal | None = None

    def two_theta(self, wavelength: float | None = None) -> np.ndarray:
        """Scattering angle 2θ in degrees.

        Parameters
        ----------
        wavelength : float, optional
            Wavelength in Å. Default: the one the pattern was calculated with.
        """
        lam = wavelength or self.settings["wavelength"]
        return np.degrees(2 * np.arcsin(np.clip(self.q * lam / (4 * np.pi), -1, 1)))

    def plot(self, ax=None, *, x: str = "q", **kwargs):
        """Plot the pattern with matplotlib.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes to draw into. Default: a new figure.
        x : {"q", "tth"}
            Plot against Q or against 2θ.
        **kwargs
            Passed on to :meth:`matplotlib.axes.Axes.plot`.

        Returns
        -------
        matplotlib.axes.Axes
        """
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots(figsize=(8, 4.5))
        if x == "q":
            xs, label = self.q, r"$Q$ (Å$^{-1}$)"
        elif x == "tth":
            xs, label = self.two_theta(), r"$2\theta$ (°)" + f",  λ = {self.settings['wavelength']} Å"
        else:
            raise ValueError("x must be 'q' or 'tth'")
        kwargs.setdefault("lw", 1.0)
        kwargs.setdefault("label", self.crystal.path.stem if self.crystal else None)
        ax.plot(xs, self.intensity, **kwargs)
        ax.set_xlabel(label)
        ax.set_ylabel("Intensity (arb. units)")
        ax.set_xlim(xs[0], xs[-1])
        return ax


def powder(crystal: Crystal, *, qmin: float = 0.5, qmax: float = 10.0, dq: float = 0.005,
           radiation: str = "xray", wavelength: float = 0.2, method: str = "debye",
           temperature: bool = True, lp_correction: str = "none") -> PowderPattern:
    """Calculate the powder pattern of a crystal with DISCUS.

    Parameters
    ----------
    crystal : Crystal
        The structure, from :func:`read_structure`.
    qmin, qmax, dq : float
        Q range and step in Å⁻¹.
    radiation : {"xray", "neutron", "electron"}
        Radiation type; selects the atomic form factors or scattering lengths.
    wavelength : float
        Wavelength in Å, used for the Lorentz-polarisation factor and for
        :meth:`PowderPattern.two_theta`.
    method : {"debye", "complete"}
        ``"debye"``: Debye equation, for any structure, including finite
        particles. ``"complete"``: Fourier transform on a reciprocal-space
        grid, for periodic crystals.
    temperature : bool
        Apply the atomic displacement parameters.
    lp_correction : {"none", "bragg", "neutron"}
        Lorentz-polarisation correction.

    Returns
    -------
    PowderPattern

    Raises
    ------
    ValueError
        For an unknown ``radiation`` or ``method``, or an invalid Q range.
    MacroError
        If the calculation fails inside DISCUS.

    Examples
    --------
    >>> p = dc.powder(cr, qmin=0.5, qmax=10, dq=0.005, radiation="neutron")  # doctest: +SKIP
    >>> p.plot()  # doctest: +SKIP
    """
    try:
        rad = _RADIATION[radiation.lower()]
        calc = _METHOD[method.lower()]
    except KeyError as exc:
        raise ValueError(f"unknown option {exc}; radiation: {sorted(set(_RADIATION))}, "
                         f"method: {sorted(_METHOD)}") from None
    if not 0 < qmin < qmax or dq <= 0:
        raise ValueError("need 0 < qmin < qmax and dq > 0")

    _ensure_loaded(crystal)
    s = session()
    out = s.workdir / "powder.dat"
    out.unlink(missing_ok=True)
    s.run(f"""
discus
  powder
    reset
    {rad}
    set axis, q
    set calc, {calc}
    set disp, off
    set delta, 0.0
    set qmin, {qmin}
    set qmax, {qmax}
    set dq, {dq}
    set profile, off
    set temp, {"use" if temperature else "ignore"}
    set wvle, {wavelength}
    set four, four
    set lpcor, {lp_correction}
    run
  exit
  output
    outf {out.name}
    value inte
    form powder, q
    run
  exit
exit
""", cwd=s.workdir)
    data = np.loadtxt(out, comments="#", ndmin=2)
    settings = dict(qmin=qmin, qmax=qmax, dq=dq, radiation=rad, wavelength=wavelength,
                    method=method, temperature=temperature, lp_correction=lp_correction)
    return PowderPattern(q=data[:, 0], intensity=data[:, 1], settings=settings, crystal=crystal)
