"""The one DISCUS suite that lives in this Python process.

For now everything runs through the f2py module ``suite_python`` (built with
``-DDIFFUSE_PYTHON=ON``): Python writes the commands to a small macro and
the suite runs it. The planned bind(C) core (migration plan, WP2) will
replace this module; the public API in ``diffuse`` is meant to stay the same.

What this layer adds on top of ``suite_python``:

* **Errors become exceptions.** A command error ends the whole macro run.
  Every macro written here sets a flag on its last line; if the flag is not
  set, ``run()`` raises ``MacroError`` with the ``***...***`` error lines.
* **Output is captured.** The suite writes to the process's stdout/stderr.
  ``run()`` sends both file descriptors to a temporary file while the macro
  runs, and keeps the text in ``Session.last_output``. Pass ``echo=True``
  (or set ``DIFFUSE_ECHO=1``) to print it as well.
* **Working directory.** The suite keeps its own working directory, so each
  macro starts with ``cd <folder>``.
"""

from __future__ import annotations

import atexit
import importlib
import itertools
import os
import re
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

# Must be set before libgfortran is loaded, so that output written while
# stdout is redirected is not flushed later, after it has been restored.
os.environ.setdefault("GFORTRAN_UNBUFFERED_PRECONNECTED", "y")

# r[] runs from 0 to MAXPAR = 500 (lib_f90/param_mod.f90); the top entry is
# reserved for the completion flag.
FLAG = "r[499]"
_ANSI = re.compile(r"\x1b\[[0-9;]*m")


class DiffuseError(RuntimeError):
    """Base class for errors raised by the diffuse package."""


class MacroError(DiffuseError):
    """A command failed inside DISCUS.

    The message holds DISCUS's ``***`` error lines; ``output`` holds all
    the text DISCUS wrote during the failed call.
    """

    def __init__(self, message: str, output: str = ""):
        super().__init__(message)
        self.output = output


def _import_suite_python():
    path = os.environ.get("DIFFUSE_PYTHON_PATH")
    if path and path not in sys.path:
        sys.path.insert(0, path)
    try:
        return importlib.import_module("suite_python").suite
    except ImportError as exc:
        raise ImportError(
            "The DISCUS Python module 'suite_python' was not found. Build "
            "DiffuseCode with -DDIFFUSE_PYTHON=ON and set "
            "DIFFUSE_PYTHON_PATH=<build>/python.") from exc


class _Capture:
    """Send file descriptors 1 and 2 to a temporary file.

    ``flush`` writes out what the Fortran runtime still buffers. Setting
    GFORTRAN_UNBUFFERED_PRECONNECTED is not enough on its own: libgfortran
    reads it once, when it is first loaded, and NumPy (through OpenBLAS) may
    have loaded it before this module.
    """

    def __init__(self, flush=None):
        self.text = ""
        self._flush = flush

    def __enter__(self):
        if self._flush is not None:
            self._flush()               # output from before the capture
        sys.stdout.flush()
        sys.stderr.flush()
        self._file = tempfile.TemporaryFile(mode="w+b")
        self._saved = [os.dup(1), os.dup(2)]
        os.dup2(self._file.fileno(), 1)
        os.dup2(self._file.fileno(), 2)
        return self

    def __exit__(self, *exc):
        if self._flush is not None:
            self._flush()
        os.dup2(self._saved[0], 1)
        os.dup2(self._saved[1], 2)
        for fd in self._saved:
            os.close(fd)
        self._file.seek(0)
        self.text = _ANSI.sub("", self._file.read().decode(errors="replace"))
        self._file.close()
        return False


class Session:
    """Wrapper around the in-process DISCUS suite. Use ``diffuse.session()``."""

    def __init__(self, echo: bool | None = None):
        self.echo = bool(os.environ.get("DIFFUSE_ECHO")) if echo is None else echo
        self.last_output = ""
        self._s = _import_suite_python()
        self._count = itertools.count()
        self._workdir = Path(tempfile.mkdtemp(prefix="diffuse_"))
        atexit.register(shutil.rmtree, self._workdir, ignore_errors=True)
        self.current_crystal = None      # the Crystal now held by DISCUS
        # flush_output() is in suite_python.f90 since the demo API; builds
        # without it rely on GFORTRAN_UNBUFFERED_PRECONNECTED alone.
        self._flush = getattr(self._s, "flush_output", None)
        with self._capture() as cap:
            self._s.initialize_suite()
        self._show(cap.text)

    def _capture(self) -> _Capture:
        return _Capture(self._flush)

    def _show(self, text: str) -> None:
        self.last_output = text
        if self.echo and text:
            print(text, end="" if text.endswith("\n") else "\n")

    @property
    def workdir(self) -> Path:
        """Scratch folder for files written by DISCUS."""
        return self._workdir

    def run(self, commands: str, cwd: str | os.PathLike | None = None) -> str:
        """Run macro text with ``cwd`` as working directory; return the output.

        Raises MacroError if a command error stopped the macro.
        """
        cwd = Path(cwd or self._workdir).resolve()
        driver = self._workdir / f"_driver_{next(self._count)}.mac"
        driver.write_text(
            f"cd {cwd}\n"
            f"{FLAG} = 0\n"
            f"{commands.strip()}\n"
            f"{FLAG} = 1\n")
        with self._capture() as cap:
            self._s.execute_macro("@" + str(driver))
            ok = self._s.get_value(FLAG) == 1.0
        self._show(cap.text)
        if not ok:
            errors = [line.strip() for line in cap.text.splitlines() if "***" in line]
            raise MacroError("DISCUS command failed:\n  " + "\n  ".join(errors or ["(no message)"]),
                             cap.text)
        return cap.text

    def value(self, expression: str) -> float:
        """Value of a DISCUS expression, e.g. ``n[1]`` or ``lat[1]``."""
        with self._capture() as cap:
            v = float(self._s.get_value(expression))
        if "***" in cap.text:
            raise MacroError(f"cannot evaluate '{expression}':\n  " + cap.text.strip(), cap.text)
        return v

    # -- KUPLOT data sets -------------------------------------------------
    def n_datasets(self) -> int:
        return int(self._s.get_data_sets())

    def dataset(self, ik: int) -> tuple[np.ndarray, np.ndarray]:
        """KUPLOT data set ``ik`` (1-based) as float32 arrays ``(x, y)``."""
        self._check_ik(ik)
        n = int(self._s.get_data_length(ik))
        x = np.empty(n, np.float32)
        y = np.empty(n, np.float32)
        self._s.get_data(ik, x, y)
        return x, y

    def dataset_2d(self, ik: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """KUPLOT grid data set ``ik`` as ``(x[nx], y[ny], z[nx, ny])``."""
        self._check_ik(ik)
        n = int(self._s.get_data_length(ik))
        x, y, z = (np.empty(n, np.float32) for _ in range(3))
        nx = np.zeros(1, np.int32)
        ny = np.zeros(1, np.int32)
        self._s.get_data_3d(ik, x, y, z, nx, ny)
        nx, ny = int(nx[0]), int(ny[0])
        if nx * ny != n:
            raise ValueError(f"data set {ik} is not a grid ({nx} x {ny} != {n})")
        return x[::ny].copy(), y[:ny].copy(), z.reshape(nx, ny)

    def _check_ik(self, ik: int) -> None:
        if not 1 <= ik <= self.n_datasets():
            raise IndexError(f"KUPLOT data set {ik} does not exist "
                             f"({self.n_datasets()} loaded)")


_SESSION: Session | None = None


def session(echo: bool | None = None) -> Session:
    """The process-wide Session (created on first use)."""
    global _SESSION
    if _SESSION is None:
        _SESSION = Session(echo=echo)
    elif echo is not None:
        _SESSION.echo = echo
    return _SESSION
