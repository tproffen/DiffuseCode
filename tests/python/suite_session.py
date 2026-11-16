"""Thin test helper around the f2py module ``suite_python``.

``suite_python`` (built with ``-DDIFFUSE_PYTHON=ON``) runs the whole
DISCUS suite inside the Python process. Its interface is small:

    initialize_suite()              once per process
    execute_macro("@/abs/path.mac") run a macro (errors are only printed)
    get_value("expr")               evaluate an expression, e.g. "n[1]"
    get_data_sets()                 number of KUPLOT data sets
    get_data_length(ik)             points in KUPLOT data set ik
    get_data(ik, x, y)              fill float32 arrays x, y in place
    get_data_3d(ik, x, y, z, nx, ny)  same for a 2-D (grid) data set

This helper adds what the tests need on top of that:

* ``run()`` writes the commands to a driver macro that first changes into a
  working folder (the suite keeps its own working directory; os.chdir does
  not change it) and sets a flag on its last line. A command error ends the
  whole macro run, so a flag that is not set means the macro failed, and
  ``run()`` raises ``MacroError``.
* ``dataset()`` / ``dataset_2d()`` return KUPLOT data sets as NumPy arrays.

The module is found through ``$DIFFUSE_PYTHON_PATH`` (the folder that holds
``suite_python*.so``, normally ``<build>/python``) or the normal import path.
"""

from __future__ import annotations

import importlib
import itertools
import os
import sys
from pathlib import Path

import numpy as np

# r[] runs from 0 to MAXPAR = 500 (lib_f90/param_mod.f90). Tests leave the
# top of the range to this helper.
FLAG = "r[499]"


class MacroError(RuntimeError):
    """A macro run through ``SuiteSession.run`` did not reach its end."""


def import_suite_python():
    """Return the ``suite`` module object of suite_python, or raise ImportError."""
    path = os.environ.get("DIFFUSE_PYTHON_PATH")
    if path and path not in sys.path:
        sys.path.insert(0, path)
    return importlib.import_module("suite_python").suite


class SuiteSession:
    """One initialised DISCUS suite in this process (there can only be one)."""

    def __init__(self, module):
        self.s = module
        self.s.initialize_suite()
        self._count = itertools.count()

    # -- running commands -------------------------------------------------
    def run(self, commands: str, cwd: Path) -> None:
        """Run ``commands`` (macro text) with ``cwd`` as working directory.

        Raises MacroError if a command error stopped the macro.
        """
        cwd = Path(cwd).resolve()
        driver = cwd / f"_pytest_driver_{next(self._count)}.mac"
        driver.write_text(
            f"cd {cwd}\n"
            f"{FLAG} = 0\n"
            f"{commands.strip()}\n"
            f"{FLAG} = 1\n")
        self.s.execute_macro("@" + str(driver))
        if self.value(FLAG) != 1.0:
            raise MacroError(
                f"macro did not finish (see the ***...*** lines in the captured "
                f"output); driver: {driver}")

    def value(self, expression: str) -> float:
        """Value of a suite expression, e.g. ``n[1]`` or ``lat[1]``."""
        return float(self.s.get_value(expression))

    # -- KUPLOT data ------------------------------------------------------
    def n_datasets(self) -> int:
        return int(self.s.get_data_sets())

    def dataset(self, ik: int) -> tuple[np.ndarray, np.ndarray]:
        """KUPLOT data set ``ik`` (1-based) as arrays ``(x, y)``."""
        self._check_ik(ik)
        n = int(self.s.get_data_length(ik))
        x = np.empty(n, np.float32)
        y = np.empty(n, np.float32)
        self.s.get_data(ik, x, y)
        return x, y

    def dataset_2d(self, ik: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """KUPLOT grid data set ``ik`` as ``(x[nx], y[ny], z[nx, ny])``."""
        self._check_ik(ik)
        n = int(self.s.get_data_length(ik))
        x = np.empty(n, np.float32)
        y = np.empty(n, np.float32)
        z = np.empty(n, np.float32)
        nx = np.zeros(1, np.int32)
        ny = np.zeros(1, np.int32)
        self.s.get_data_3d(ik, x, y, z, nx, ny)
        nx, ny = int(nx[0]), int(ny[0])
        if nx * ny != n:
            raise ValueError(f"data set {ik} is not a grid ({nx} x {ny} != {n})")
        return x[::ny].copy(), y[:ny].copy(), z.reshape(nx, ny)

    def _check_ik(self, ik: int) -> None:
        if not 1 <= ik <= self.n_datasets():
            raise IndexError(f"KUPLOT data set {ik} does not exist "
                             f"({self.n_datasets()} loaded)")
