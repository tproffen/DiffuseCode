"""The one DISCUS suite that lives in this Python process.

The compiled core (``_core.py``) is the C library ``libdiffuse_core``,
loaded with ctypes (built with ``-DDIFFUSE_PYTHON=ON``, or by ``pip install``). Python hands the commands to the suite as a
macro held in memory (no file), and the suite runs it.

What this layer adds on top of the core:

* **Errors become exceptions.** A command error ends the whole macro run.
  Every macro written here sets a flag on its last line; if the flag is not
  set, ``run()`` raises ``MacroError``. DISCUS runs in library mode
  (``dc_init``): an error never ends the process, and the first error is
  recorded, so ``MacroError`` carries its number, type and message.
* **Output is captured.** The suite writes to the process's stdout/stderr
  (unbuffered, so line by line). ``run()`` sends both file descriptors to a
  temporary file while the macro runs, and keeps the text in
  ``Session.last_output``. With ``echo=True`` (or ``DIFFUSE_ECHO=1``) the
  lines are also shown as DISCUS writes them, in the terminal or in the
  notebook cell: the call into the core releases the GIL, and a thread follows the
  file. In library mode DISCUS keeps all of its output there: ``set prompt``
  cannot move it to a file, and errors write all their lines to it.
* **Working directory.** The suite keeps its own working directory, so each
  macro starts with ``cd <folder>``.
"""

from __future__ import annotations

import atexit
import os
import re
import shutil
import sys
import tempfile
import threading
from pathlib import Path

import numpy as np

from . import _core

# Must be set before libgfortran is loaded, so that output written while
# stdout is redirected is not flushed later, after it has been restored.
os.environ.setdefault("GFORTRAN_UNBUFFERED_PRECONNECTED", "y")

# r[] runs from 0 to MAXPAR = 500 (lib_f90/param_mod.f90); the top entry is
# reserved for the completion flag.
FLAG = "r[499]"
# Longest command line DISCUS takes (PREC_STRING in lib_f90/precision_mod.f90)
MAX_LINE = 1024
_ANSI = re.compile(r"\x1b\[[0-9;]*m")


class DiffuseError(RuntimeError):
    """Base class for errors raised by the diffuse package."""


class MacroError(DiffuseError):
    """A command failed inside DISCUS.

    Attributes
    ----------
    code : int
        DISCUS error number (negative), or 0 if DISCUS recorded none.
    category : str
        DISCUS error type: ``"COMM"`` (command language), ``"IO"``,
        ``"FORT"`` (expressions), ``"MAC"``, ``"MATH"``, ``"APPL"`` (the
        program itself) and others.
    text : str
        DISCUS's standard message for this error number.
    details : list of str
        Further lines DISCUS gave for this error (often the file or value
        concerned).
    output : str
        All the text DISCUS wrote during the failed call.
    """

    def __init__(self, message: str, output: str = "", code: int = 0,
                 category: str = "", text: str = "", details=()):
        super().__init__(message)
        self.output = output
        self.code = code
        self.category = category
        self.text = text
        self.details = list(details)


class _Capture:
    """Send file descriptors 1 and 2 to a temporary file.

    ``flush`` writes out what the Fortran runtime still buffers. Setting
    GFORTRAN_UNBUFFERED_PRECONNECTED is not enough on its own: libgfortran
    reads it once, when it is first loaded, and NumPy (through OpenBLAS) may
    have loaded it before this module.

    ``stdin``: a file of empty lines that becomes file descriptor 0 for the
    call, so that DISCUS never waits for the keyboard: a read gets an empty
    line, as if Return was pressed (a pause continues, a question keeps its
    default). Not /dev/null: after an end of file, gfortran stops the whole
    process at the next read of that unit.

    With ``live=True`` a thread shows complete lines while they are written
    (on the real terminal, or through ``sys.stdout`` in Jupyter, where it is
    not the process's file descriptor 1); ``shown`` is the text shown so far.
    """

    POLL = 0.2                          # seconds between looks at the file

    def __init__(self, flush=None, live: bool = False, stdin=None):
        self.text = ""
        self.shown = ""
        self._flush = flush
        self._live = live
        self._stdin = stdin

    def __enter__(self):
        if self._flush is not None:
            self._flush()               # output from before the capture
        sys.stdout.flush()
        sys.stderr.flush()
        self._file = tempfile.TemporaryFile(mode="w+b")
        self._saved = [os.dup(1), os.dup(2)]
        os.dup2(self._file.fileno(), 1)
        os.dup2(self._file.fileno(), 2)
        self._saved_in = False                # False: not replaced; None: fd 0 was closed
        if self._stdin is not None:
            try:
                self._saved_in = os.dup(0)
            except OSError:
                self._saved_in = None
            fd = os.open(self._stdin, os.O_RDONLY)
            os.dup2(fd, 0)
            os.close(fd)
        self._pos = 0
        self._thread = None
        if self._live:
            self._stop = threading.Event()
            self._thread = threading.Thread(target=self._follow, daemon=True)
            self._thread.start()
        return self

    def __exit__(self, *exc):
        if self._flush is not None:
            self._flush()
        if self._thread is not None:
            self._stop.set()
            self._thread.join()
            self._show_new()            # what came in since the last look
        os.dup2(self._saved[0], 1)
        os.dup2(self._saved[1], 2)
        for fd in self._saved:
            os.close(fd)
        if self._saved_in is None:
            os.close(0)
        elif self._saved_in is not False:
            os.dup2(self._saved_in, 0)
            os.close(self._saved_in)
        self._file.seek(0)
        self.text = _ANSI.sub("", self._file.read().decode(errors="replace"))
        self._file.close()
        return False

    # -- live output -------------------------------------------------------
    def _follow(self):
        while not self._stop.wait(self.POLL):
            self._show_new()

    def _show_new(self):
        """Show the complete lines written since the last call."""
        size = os.fstat(self._file.fileno()).st_size
        if size <= self._pos:
            return
        data = os.pread(self._file.fileno(), size - self._pos, self._pos)
        end = data.rfind(b"\n") + 1
        if not end:
            return
        self._pos += end
        text = _ANSI.sub("", data[:end].decode(errors="replace"))
        self.shown += text
        if sys.stdout is sys.__stdout__:    # a terminal or pipe: fd 1 is the capture now
            os.write(self._saved[0], text.encode())
        else:                               # Jupyter and similar: their own stream
            sys.stdout.write(text)
            sys.stdout.flush()


class Session:
    """Wrapper around the in-process DISCUS suite. Use ``diffuse.session()``."""

    def __init__(self, echo: bool | None = None):
        self.echo = bool(os.environ.get("DIFFUSE_ECHO")) if echo is None else echo
        self.last_output = ""
        self._workdir = Path(tempfile.mkdtemp(prefix="diffuse_"))
        atexit.register(shutil.rmtree, self._workdir, ignore_errors=True)
        self._c = _core.load()
        # Standard input for DISCUS calls (see _Capture): empty lines
        self._stdin = self._workdir / "_stdin"
        self._stdin.write_text("\n" * 20000)
        self.current_crystal = None      # the Crystal now held by DISCUS
        with self._capture() as cap:        # library mode: errors as records, no STOP
            self._c.init()
        self._show(cap.text)

    @property
    def backend(self) -> str:
        """How DISCUS is reached: ``"ctypes"``, the C library libdiffuse_core
        (earlier builds also had ``"f2py"``, the module suite_python)."""
        return self._c.name

    def _capture(self, live: bool = False) -> _Capture:
        return _Capture(self._c.flush, live=live and self.echo, stdin=self._stdin)

    def _show(self, text: str, shown: str = "") -> None:
        """Keep ``text`` as the last output; with echo, print what was not
        shown live already."""
        self.last_output = text
        rest = text[len(shown):] if text.startswith(shown) else text
        if self.echo and rest:
            print(rest, end="" if rest.endswith("\n") else "\n")

    @property
    def workdir(self) -> Path:
        """Scratch folder for files written by DISCUS."""
        return self._workdir

    def run(self, commands: str, cwd: str | os.PathLike | None = None) -> str:
        """Run macro text with ``cwd`` as working directory; return the output.

        Raises MacroError if a command error stopped the macro.
        """
        cwd = Path(cwd or self._workdir).resolve()
        text = f"cd {cwd}\n{FLAG} = 0\n{commands.strip()}\n{FLAG} = 1\n"
        too_long = [i for i, line in enumerate(text.splitlines(), -1)
                    if len(line.encode()) > MAX_LINE]
        if too_long:
            raise ValueError(f"command line {too_long[0]} is longer than DISCUS's limit "
                             f"of {MAX_LINE} characters")
        # DISCUS's "cd" changes the working directory of the whole process;
        # give the caller theirs back afterwards.
        caller_cwd = os.getcwd()
        self._clear_error()
        try:
            with self._capture(live=True) as cap:
                self._c.run(text)
                ok = self._c.value(FLAG) == 1.0
        finally:
            os.chdir(caller_cwd)
        self._show(cap.text, cap.shown)
        if not ok:
            raise self._error("DISCUS command failed", cap.text)
        return cap.text

    # -- errors -------------------------------------------------------------
    def _clear_error(self) -> None:
        self._c.clear_error()

    def last_error(self) -> tuple[int, str, str, list[str]]:
        """The first error DISCUS recorded in the last call.

        Returns ``(code, category, text, details)``; ``code`` is 0 when
        there was no error.
        """
        num, cls, text, msgs = self._c.last_error()
        details = [msgs[i:i + 80].strip() for i in range(0, len(msgs), 80)]
        return int(num), cls.strip(), text.strip(), [d for d in details if d]

    def _error(self, what: str, output: str) -> MacroError:
        code, category, text, details = self.last_error()
        lines = [line.strip() for line in output.splitlines() if "***" in line]
        if code:
            # The recorded error says it all; the other *** lines are the
            # macro and menu levels reporting that it happened.
            head = f"{what}: {text} (DISCUS error {code}, {category})"
            body = details
        else:
            head, body = what + ":", lines or [
                "DISCUS stopped before the end of the commands without an error "
                "(for example 'exit' at the top level, which ends the macro)"]
        return MacroError("\n  ".join([head, *body]), output, code, category, text, details)

    def value(self, expression: str) -> float:
        """Value of a DISCUS expression, e.g. ``n[1]`` or ``lat[1]``."""
        self._clear_error()
        with self._capture() as cap:
            v = self._c.value(expression)
        if self.last_error()[0] < 0 or "***" in cap.text:
            raise self._error(f"cannot evaluate '{expression}'", cap.text)
        return v

    # -- state ---------------------------------------------------------------
    def reset(self) -> None:
        """Put DISCUS, KUPLOT, DIFFEV and REFINE back into their start-up state.

        Same as the suite command ``reset``. Variables, the random number
        generator and settings made with ``set`` (``set parallel``,
        ``set error`` ...) are kept. ``current_crystal`` is cleared, so the
        next calculation on a Crystal reads it in again.
        """
        with self._capture() as cap:
            self._c.reset()
        self._show(cap.text)
        self.current_crystal = None

    def seed(self, value: int) -> None:
        """Seed the random number generator, as the DISCUS command ``seed``.

        The same seed gives the same sequence of random numbers, so
        calculations with random elements (``lots``, MMC, ...) repeat.
        ``value=0`` seeds from the clock. With more than one OpenMP thread
        the order in which threads draw numbers can still vary.
        """
        self._c.seed(int(value))

    def random_state(self) -> np.ndarray:
        """The random number generator's complete state, for :meth:`set_random_state`.

        Includes the Gaussian number that DISCUS may keep waiting (Gaussian
        numbers come in pairs), so a restore repeats the sequence exactly.
        """
        return self._c.random_state()

    def set_random_state(self, state) -> None:
        """Put the generator back into a state from :meth:`random_state`."""
        state = np.asarray(state, dtype=np.int32)
        if state.shape != (self._c.random_size(),):
            raise ValueError(f"state must have {self._c.random_size()} integers, "
                             f"not shape {state.shape}")
        self._c.set_random_state(state)

    # -- KUPLOT data sets -------------------------------------------------
    def n_datasets(self) -> int:
        return self._c.data_sets()

    def dataset(self, ik: int) -> tuple[np.ndarray, np.ndarray]:
        """KUPLOT data set ``ik`` (1-based) as float32 arrays ``(x, y)``."""
        self._check_ik(ik)
        return self._c.data(ik)

    def dataset_2d(self, ik: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """KUPLOT grid data set ``ik`` as ``(x[nx], y[ny], z[nx, ny])``."""
        self._check_ik(ik)
        n = self._c.data_length(ik)
        try:
            x, y, z, nx, ny = self._c.data_2d(ik)
        except ValueError:
            nx = ny = 0
        if nx * ny != n:
            raise ValueError(f"data set {ik} is not a grid ({nx} x {ny} != {n})")
        return x[::ny].copy(), y[:ny].copy(), z.reshape(nx, ny)

    def _check_ik(self, ik: int) -> None:
        if not 1 <= ik <= self.n_datasets():
            raise IndexError(f"KUPLOT data set {ik} does not exist "
                             f"({self.n_datasets()} loaded)")


_SESSION: Session | None = None


def session(echo: bool | None = None) -> Session:
    """The process-wide Session (created on first use).

    Parameters
    ----------
    echo : bool, optional
        Show DISCUS's output as it is written (in a terminal or a Jupyter
        cell), for example to follow a long calculation. Default: off, or
        on if the environment variable ``DIFFUSE_ECHO`` is set. The output
        of the last call is always in ``session().last_output``.
    """
    global _SESSION
    if _SESSION is None:
        _SESSION = Session(echo=echo)
    elif echo is not None:
        _SESSION.echo = echo
    return _SESSION
