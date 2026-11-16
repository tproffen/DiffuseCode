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
    """A calculated powder pattern on an evenly spaced Q grid."""

    q: np.ndarray                   # Å^-1
    intensity: np.ndarray
    settings: dict = field(default_factory=dict)
    crystal: Crystal | None = None

    def two_theta(self, wavelength: float | None = None) -> np.ndarray:
        """2θ in degrees for ``wavelength`` (default: the one used)."""
        lam = wavelength or self.settings["wavelength"]
        return np.degrees(2 * np.arcsin(np.clip(self.q * lam / (4 * np.pi), -1, 1)))

    def plot(self, ax=None, *, x: str = "q", **kwargs):
        """Plot with matplotlib; ``x`` is ``"q"`` or ``"tth"``. Returns the Axes."""
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
    """Calculate the powder pattern of ``crystal`` with DISCUS.

    Parameters
    ----------
    qmin, qmax, dq : Q range and step in Å^-1.
    radiation : "xray", "neutron" or "electron".
    wavelength : in Å (used for the Q limit, Lorentz-polarisation and 2θ).
    method : "debye" (Debye equation, any structure) or "complete"
        (Fourier on a reciprocal-space grid, for periodic crystals).
    temperature : apply the atomic displacement parameters.
    lp_correction : "none", "bragg" or "neutron" Lorentz-polarisation factor.
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
