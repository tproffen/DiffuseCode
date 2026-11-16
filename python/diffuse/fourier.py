"""Single-crystal diffraction: intensities on a line or plane in reciprocal space."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ._session import session
from .crystal import Crystal, _ensure_loaded

_RADIATION = {"xray": "xray", "x-ray": "xray", "neutron": "neutron", "electron": "electron"}
_HKL = "hkl"


@dataclass
class FourierMap:
    """Intensities on a line or plane in reciprocal space. See :func:`fourier`.

    The points are ``ll + i/(na-1) * (lr - ll) + j/(no-1) * (ul - ll)`` for
    ``i < na`` along the abscissa and ``j < no`` along the ordinate.

    Attributes
    ----------
    intensity : numpy.ndarray
        Shape ``(na, no)`` for a plane, ``(na,)`` for a line; index 0 runs
        along the abscissa (``ll`` to ``lr``).
    ll, lr, ul : numpy.ndarray
        Corners in reciprocal-lattice units (h, k, l): lower left, lower
        right, upper left.
    x, y : numpy.ndarray
        Axis values used for plotting: the Miller index that changes most
        along each axis (named in ``x_label``, ``y_label``). ``y`` is empty
        for a line.
    settings : dict
        The arguments the map was calculated with.
    crystal : Crystal or None
        The structure it was calculated for.
    """

    intensity: np.ndarray
    ll: np.ndarray
    lr: np.ndarray
    ul: np.ndarray
    x: np.ndarray
    y: np.ndarray
    x_label: str
    y_label: str
    settings: dict = field(default_factory=dict)
    crystal: Crystal | None = None

    @property
    def is_line(self) -> bool:
        """True for a 1-D line, False for a 2-D plane."""
        return self.intensity.ndim == 1

    def hkl(self) -> np.ndarray:
        """Reciprocal-space coordinates (h, k, l) of every point.

        Returns
        -------
        numpy.ndarray
            Shape ``(na, 3)`` for a line, ``(na, no, 3)`` for a plane.
        """
        na = self.intensity.shape[0]
        u = np.linspace(0, 1, na)
        if self.is_line:
            return self.ll + u[:, None] * (self.lr - self.ll)
        no = self.intensity.shape[1]
        v = np.linspace(0, 1, no)
        return (self.ll + u[:, None, None] * (self.lr - self.ll)
                + v[None, :, None] * (self.ul - self.ll))

    def plot(self, ax=None, *, log: bool = False, dynamic_range: float = 1e5,
             cmap: str = "viridis", **kwargs):
        """Plot the map (image) or line with matplotlib.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes to draw into. Default: a new figure.
        log : bool
            Logarithmic intensity scale (zero intensities are masked).
        dynamic_range : float
            With ``log``, show intensities down to ``max / dynamic_range``.
            Finite crystals have fringes many decades below the peaks;
            without a limit they would swamp the map.
        cmap : str
            Colour map for a plane.
        **kwargs
            Passed on to :meth:`~matplotlib.axes.Axes.imshow` (plane) or
            :meth:`~matplotlib.axes.Axes.plot` (line).

        Returns
        -------
        matplotlib.axes.Axes
        """
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm

        if ax is None:
            _, ax = plt.subplots(figsize=(6, 5) if not self.is_line else (8, 4.5))
        if self.is_line:
            kwargs.setdefault("lw", 1.0)
            ax.plot(self.x, self.intensity, **kwargs)
            ax.set_xlabel(self.x_label)
            ax.set_ylabel("Intensity (arb. units)")
            ax.set_xlim(self.x[0], self.x[-1])
            if log:
                ax.set_yscale("log")
                top = float(self.intensity.max())
                ax.set_ylim(top / dynamic_range, top * 2)
            return ax
        z = self.intensity.T                       # rows = ordinate for imshow
        if log:
            z = np.ma.masked_less_equal(z, 0)
            top = float(z.max())
            kwargs.setdefault("norm", LogNorm(vmin=top / dynamic_range, vmax=top))
        im = ax.imshow(z, origin="lower", cmap=cmap, aspect="auto",
                       extent=(self.x[0], self.x[-1], self.y[0], self.y[-1]), **kwargs)
        ax.set_xlabel(self.x_label)
        ax.set_ylabel(self.y_label)
        ax.figure.colorbar(im, ax=ax, label="Intensity")
        return ax


def _axis_index(vector: np.ndarray, taken: int | None = None) -> int:
    """Miller index (0=h, 1=k, 2=l) that changes most along ``vector``."""
    order = np.argsort(-np.abs(vector), kind="stable")
    for i in order:
        if i != taken:
            return int(i)
    return int(order[0])


def fourier(crystal: Crystal, *, ll, lr, ul=None, npoints=(201, 201),
            radiation: str = "xray", wavelength: float | None = None,
            temperature: bool = True, subtract_average: float = 0.0) -> FourierMap:
    """Calculate single-crystal intensities on a line or plane with DISCUS.

    Parameters
    ----------
    crystal : Crystal
        The structure, from :func:`read_structure`.
    ll, lr : sequence of 3 float
        Lower-left and lower-right corner (h, k, l); the abscissa runs from
        ``ll`` to ``lr``.
    ul : sequence of 3 float, optional
        Upper-left corner; the ordinate runs from ``ll`` to ``ul``. Leave it
        out for a line.
    npoints : int or (int, int)
        Number of points along the abscissa and the ordinate (a plain int
        for a line, or the same along both axes of a plane).
    radiation : {"xray", "neutron", "electron"}
        Radiation type; selects the atomic form factors or scattering lengths.
    wavelength : float, optional
        Wavelength in Å. Default: DISCUS's current setting.
    temperature : bool
        Apply the atomic displacement parameters.
    subtract_average : float
        Percentage (0–100) of the average structure's scattering to subtract.
        100 removes the Bragg peaks of the average structure and leaves the
        diffuse scattering of a disordered crystal.

    Returns
    -------
    FourierMap

    Raises
    ------
    ValueError
        For an unknown ``radiation``, bad corners or point counts.
    MacroError
        If the calculation fails inside DISCUS.

    Examples
    --------
    The hk0 plane from -3 to 3 in h and k:

    >>> f = dc.fourier(cr, ll=(-3, -3, 0), lr=(3, -3, 0), ul=(-3, 3, 0),
    ...                npoints=(301, 301))  # doctest: +SKIP
    >>> f.plot(log=True)  # doctest: +SKIP
    """
    try:
        rad = _RADIATION[radiation.lower()]
    except KeyError:
        raise ValueError(f"unknown radiation '{radiation}'; use {sorted(set(_RADIATION))}") from None
    ll = np.asarray(ll, dtype=float)
    lr = np.asarray(lr, dtype=float)
    line = ul is None
    ul = ll.copy() if line else np.asarray(ul, dtype=float)
    if ll.shape != (3,) or lr.shape != (3,) or ul.shape != (3,):
        raise ValueError("ll, lr and ul must each be three numbers (h, k, l)")
    if np.allclose(lr, ll) or (not line and np.allclose(ul, ll)):
        raise ValueError("lr (and ul for a plane) must differ from ll")
    if line:
        na, no = (int(npoints), 1) if np.ndim(npoints) == 0 else (int(npoints[0]), 1)
    else:
        na, no = (int(npoints),) * 2 if np.ndim(npoints) == 0 else map(int, npoints)
    if na < 2 or (not line and no < 2):
        raise ValueError("need at least 2 points along each axis")
    if not 0.0 <= subtract_average <= 100.0:
        raise ValueError("subtract_average is a percentage, 0 to 100")

    ia = _axis_index(lr - ll)
    io = _axis_index(ul - ll, taken=ia) if not line else (ia + 1) % 3

    def vec(v):
        return ", ".join(f"{c:.6g}" for c in v)

    _ensure_loaded(crystal)
    s = session()
    out = s.workdir / "fourier.inte"
    out.unlink(missing_ok=True)
    wvle = f"\n    wvle {wavelength}" if wavelength is not None else ""
    s.run(f"""
discus
  fourier
    reset
    {rad}
    ll {vec(ll)}
    lr {vec(lr)}
    ul {vec(ul)}
    na {na}
    no {no}
    abs {_HKL[ia]}
    ord {_HKL[io]}
    temp {"use" if temperature else "ign"}
    disp off{wvle}
    set aver, {subtract_average}
    lots off
    run
  exit
  output
    outf {out.name}
    value inte
    form stan
    run
  exit
exit
kuplot
  rese
  load {"xy" if line else "ni"}, {out.name}
exit
""", cwd=s.workdir)
    if line:
        x, inte = s.dataset(1)
        y = np.empty(0)
        intensity = inte.astype(float)
    else:
        x, y, z = s.dataset_2d(1)
        intensity = z.astype(float)
    settings = dict(ll=tuple(ll), lr=tuple(lr), ul=None if line else tuple(ul),
                    npoints=(na,) if line else (na, no), radiation=rad,
                    wavelength=wavelength, temperature=temperature,
                    subtract_average=subtract_average)
    return FourierMap(intensity=intensity, ll=ll, lr=lr, ul=ul,
                      x=np.asarray(x, dtype=float), y=np.asarray(y, dtype=float),
                      x_label=_HKL[ia], y_label="" if line else _HKL[io],
                      settings=settings, crystal=crystal)
