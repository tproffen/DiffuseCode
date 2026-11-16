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
from .crystal import Crystal, read_structure
from .powder import PowderPattern, powder

__all__ = [
    "Crystal", "DiffuseError", "MacroError", "PowderPattern", "Session",
    "powder", "read_structure", "run_macro", "session",
]


def run_macro(commands: str, cwd=None) -> str:
    """Run DISCUS macro text (escape hatch for commands not wrapped yet)."""
    return session().run(commands, cwd)
