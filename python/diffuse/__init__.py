"""DiffuseCode from Python (early demo API).

    import diffuse as dc

    cr = dc.read_structure("anatase.cif", supercell=(6, 6, 3))
    p = dc.powder(cr, qmin=0.5, qmax=10, dq=0.005, radiation="xray")
    p.plot()

This first version drives DISCUS through the f2py module ``suite_python``
(build with ``-DDIFFUSE_PYTHON=ON`` and set ``DIFFUSE_PYTHON_PATH`` to
``<build>/python``). DISCUS holds one structure at a time, in this process.
See migration/DiffuseCode_Python_Migration_Plan.md for where this is going.
"""

from ._session import DiffuseError, MacroError, Session, session
from .crystal import Crystal, StructureView, read_structure
from .fourier import FourierMap, fourier
from .pdf import PDFPattern, pdf
from .powder import PowderPattern, powder

__all__ = [
    "Crystal", "DiffuseError", "FourierMap", "MacroError", "PDFPattern",
    "PowderPattern", "Session", "StructureView",
    "fourier", "pdf", "powder", "read_structure", "reset", "run_macro", "seed",
    "session",
]


def run_macro(commands: str, cwd=None) -> str:
    """Run DISCUS macro text (escape hatch for commands not wrapped yet)."""
    return session().run(commands, cwd)


def seed(value: int) -> None:
    """Seed DISCUS's random number generator.

    Calculations with random elements (``lots`` in Fourier, Monte Carlo,
    random displacements) then repeat exactly. ``value=0`` seeds from the
    clock.

    Parameters
    ----------
    value : int
        The seed, as for the DISCUS command ``seed``.

    Examples
    --------
    >>> import diffuse as dc
    >>> dc.seed(1234)
    """
    session().seed(value)


def reset() -> None:
    """Put DISCUS back into its start-up state (the suite command ``reset``).

    Removes the structure and all menu settings. Variables, the random
    number generator and settings made with ``set`` (for example
    ``set parallel``) are kept. A :class:`Crystal` made earlier is read in
    again when it is used next.
    """
    session().reset()
