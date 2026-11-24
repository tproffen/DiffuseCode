"""WP2 core: library mode, errors as exceptions, reset and random state.

The "done when" test of WP2 in the migration plan: a deliberately bad input
raises a Python exception, the process survives, nothing reaches the
terminal, and the next call works.

    DIFFUSE_PYTHON_PATH=<build>/python pytest tests/python/test_wp2_core.py
"""

from __future__ import annotations

import os
import sys
import time
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[2]
NI = REPO / "python" / "examples" / "nickel.cell"


@pytest.fixture
def dc(suite):
    import diffuse
    return diffuse


# ------------------------------------------------------------------ errors
def test_bad_command_raises_with_discus_error(dc, capfd):
    capfd.readouterr()
    with pytest.raises(dc.MacroError) as info:
        dc.run_macro("discus\n  no_such_command 1\nexit")
    err = info.value
    assert (err.code, err.category, err.text) == (-8, "COMM", "Unknown command")
    assert "DISCUS error -8" in str(err)
    out, errout = capfd.readouterr()
    assert out == "" and errout == ""                  # nothing on the terminal
    assert dc.session().value("2 * 21") == 42          # and the next call works


def test_missing_file_error_details(dc):
    with pytest.raises(dc.MacroError) as info:
        dc.run_macro("discus\n  read\n    stru /no/such/file.stru\n  exit\nexit")
    err = info.value
    assert err.code < 0 and err.category == "I/O"
    assert "/no/such/file.stru" in err.details


def test_set_error_exit_does_not_end_the_process(dc):
    """'set error,exit' made the first error STOP the program (and Python
    with it). In library mode it is an error like any other."""
    try:
        with pytest.raises(dc.MacroError) as info:
            dc.run_macro("set error,exit\ndiscus\n  no_such_command\nexit")
        assert info.value.code == -8
        assert dc.session().value("1 + 1") == 2
    finally:
        dc.run_macro("set error,continue")


def test_first_error_is_reported(dc):
    """Later 'returned with error' lines of the menu levels do not replace
    the error that started it."""
    with pytest.raises(dc.MacroError) as info:
        dc.run_macro("discus\n  read\n    no_such_read_command\n  exit\nexit")
    assert info.value.text != "Sub menu returned with error"
    assert dc.session().last_error()[0] == info.value.code


def test_value_error(dc):
    with pytest.raises(dc.MacroError) as info:
        dc.session().value("no_such_variable[3]")
    assert (info.value.code, info.value.category) == (-2, "FORT")
    assert dc.session().value("3") == 3


def test_successful_call_clears_the_error(dc):
    with pytest.raises(dc.MacroError):
        dc.run_macro("discus\n  no_such_command\nexit")
    dc.run_macro("r[1] = 5")
    assert dc.session().last_error() == (0, "", "", [])


# ------------------------------------------------------------------- reset
def test_reset_removes_structure_and_reloads_on_use(dc):
    cr = dc.read_structure(NI, supercell=2)
    first = dc.powder(cr, qmax=5.0)
    dc.reset()
    dc.run_macro("discus\n  r[1] = n[1]\nexit")
    assert dc.session().value("r[1]") == 0             # no atoms after reset
    assert dc.session().current_crystal is None
    again = dc.powder(cr, qmax=5.0)                    # read in again
    np.testing.assert_allclose(again.intensity, first.intensity, rtol=1e-10)


def test_reset_keeps_variables(dc):
    dc.run_macro("r[7] = 3.5")
    dc.reset()
    assert dc.session().value("r[7]") == 3.5


# ----------------------------------------------------------- random numbers
def draws(s, n=5):
    return [s.value("ran(0)") for _ in range(n)]


def test_seed_repeats_and_differs(dc):
    s = dc.session()
    dc.seed(1234)
    a = draws(s)
    dc.seed(1234)
    assert draws(s) == a
    dc.seed(4321)
    assert draws(s) != a


def test_seed_matches_seed_command(dc):
    s = dc.session()
    dc.seed(77)
    a = draws(s)
    dc.run_macro("seed 77")
    assert draws(s) == a


def test_random_state_round_trip(dc):
    s = dc.session()
    dc.seed(99)
    draws(s, 3)
    state = s.random_state()
    assert state.dtype.kind == "i" and state.ndim == 1
    a = draws(s) + [s.value("gran(1.0, s)") for _ in range(3)]
    s.set_random_state(state)
    b = draws(s) + [s.value("gran(1.0, s)") for _ in range(3)]
    assert a == b


def test_random_state_bad_shape(dc):
    with pytest.raises(ValueError):
        dc.session().set_random_state([1, 2])


# ------------------------------------------------------------------ output
class TimedStream:
    """Stands in for sys.stdout and notes when each text arrives."""

    def __init__(self):
        self.parts = []

    def write(self, text):
        self.parts.append((time.monotonic(), text))
        return len(text)

    def flush(self):
        pass

    def lines_with(self, word):
        return [(t, line) for t, text in self.parts for line in text.splitlines()
                if word in line]


def test_echo_shows_output_while_the_macro_runs(dc, monkeypatch):
    """With echo, lines appear as DISCUS writes them (the ctypes call releases
    the GIL and a thread follows the output), and each line only once."""
    s = dc.session()
    stream = TimedStream()
    monkeypatch.setattr(sys, "stdout", stream)
    monkeypatch.setattr(s, "echo", True)
    start = time.monotonic()
    dc.run_macro('do i[1] = 1, 3\n  echo "live step %d", i[1]\n  system "sleep 0.6"\nenddo')
    end = time.monotonic()
    steps = [(t, line) for t, line in stream.lines_with("live step") if "echo" not in line]
    assert [line.strip() for _, line in steps] == ["live step 1", "live step 2", "live step 3"]
    times = [t for t, _ in steps]
    assert times[0] < end - 0.8                    # step 1 shown well before the end
    assert times[1] - times[0] > 0.3 and times[2] - times[1] > 0.3
    assert s.last_output.count("live step 2") == 1     # kept once in the output


def test_without_echo_nothing_is_shown(dc, capfd):
    capfd.readouterr()
    dc.run_macro('echo "quiet"')
    assert capfd.readouterr() == ("", "")
    assert "quiet" in dc.session().last_output


def test_set_prompt_cannot_move_the_output(dc):
    """In library mode 'set prompt ..., file' (output to <program>.log)
    is ignored, so the output is still captured."""
    try:
        out = dc.run_macro('set prompt,off,file\necho "still here"')
        assert "still here" in out
        assert not list(dc.session().workdir.glob("*.log"))
    finally:
        dc.run_macro("set prompt,on")


def test_error_output_has_all_lines(dc):
    """The extra lines of an error go to the same output as the error."""
    with pytest.raises(dc.MacroError) as info:
        dc.run_macro("discus\n  read\n    stru /no/such/file.stru\n  exit\nexit")
    assert "/no/such/file.stru" in info.value.output


# ---------------------------------------------------------------- keyboard
@pytest.fixture
def typed_stdin():
    """Put text on file descriptor 0, as if typed, for the duration of a test."""
    r, w = os.pipe()
    os.write(w, b"typed by someone\n")
    os.close(w)
    saved = os.dup(0)
    os.dup2(r, 0)
    os.close(r)
    yield
    os.dup2(saved, 0)
    os.close(saved)


def test_discus_never_reads_the_keyboard(dc, typed_stdin):
    """During a DISCUS call, standard input is a supply of empty lines: what
    is waiting on the real standard input is not read (and nothing hangs)."""
    out = dc.run_macro('system "read x; echo got-[$x]"')
    assert "got-[]" in out
    assert os.read(0, 100) == b"typed by someone\n"        # still there


@pytest.mark.parametrize("macro", ["wait", "wait return"])
def test_wait_continues(dc, macro):
    start = time.monotonic()
    dc.run_macro(macro)
    assert time.monotonic() - start < 5


def test_wait_input_needs_a_terminal(dc):
    with pytest.raises(dc.MacroError) as info:
        dc.run_macro("wait input")
    assert (info.value.code, info.value.category) == (-21, "COMM")
    assert "terminal" in info.value.text
    assert dc.session().value("1 + 1") == 2


def test_top_level_exit_ends_only_the_macro(dc):
    """'exit' at the top ends the macro: an error for the call, but no
    program exit (no clean-up, no update questions); the session goes on."""
    with pytest.raises(dc.MacroError) as info:
        dc.run_macro('echo "before"\nexit\necho "after"')
    assert "stopped before the end" in str(info.value)
    assert "before" in info.value.output and " after" not in info.value.output
    assert dc.session().value("2 + 2") == 4


# ------------------------------------------------------- commands in memory
def test_no_macro_file_is_written(dc):
    s = dc.session()
    dc.run_macro('r[1] = 3')
    assert not [p.name for p in s.workdir.iterdir() if p.suffix == ".mac"]


def test_macro_files_still_run_from_commands(dc, tmp_path):
    (tmp_path / "add.mac").write_text("r[2] = $1 + $2\n")
    dc.run_macro("@add.mac 20, 22", cwd=tmp_path)
    assert dc.session().value("r[2]") == 42


def test_loops_and_blocks_in_memory(dc):
    out = dc.run_macro('r[3] = 0\ndo i[1] = 1, 5\n  r[3] = r[3] + i[1]\nenddo\n'
                       'if (r[3] == 15) then\n  echo "sum ok"\nendif')
    assert dc.session().value("r[3]") == 15 and "sum ok" in out


def test_line_too_long(dc):
    with pytest.raises(ValueError, match="longer than"):
        dc.run_macro("echo \"" + "x" * 1100 + "\"")


def test_random_state_restores_mid_gaussian_pair(dc):
    """gasdev makes Gaussian numbers in pairs; the state includes the one
    that is waiting, so a restore between the two repeats exactly."""
    s = dc.session()
    dc.seed(5)
    s.value("gran(1.0, s)")                      # first of a pair: the second waits
    state = s.random_state()
    a = [s.value("gran(1.0, s)") for _ in range(4)]
    s.set_random_state(state)
    assert [s.value("gran(1.0, s)") for _ in range(4)] == a


# ------------------------------------------------------------------ reset
RESET_SCRIPT = """
import json, sys
import numpy as np
import diffuse as dc
ex = sys.argv[1]
dc.run_macro("set parallel, useomp:serial")      # same summation order as the test
layer = dc.read_structure(ex + "/cuau_layer.stru")
f = dc.fourier(layer, ll=(0, 0, 0), lr=(2, 0, 0), ul=(0, 2, 0), npoints=(41, 41))
ni = dc.read_structure(ex + "/nickel.cell", supercell=3)
p = dc.powder(ni, qmax=6.0)
np.savez(sys.argv[2], fourier=f.intensity, powder=p.intensity)
"""


def test_reset_after_monte_carlo_matches_a_fresh_process(dc, tmp_path):
    """After an MMC run (tests/04.SRO/SRO_3, cut before its Fourier) and
    reset(), a Fourier map and a powder pattern are identical to those of a
    fresh Python process: reset leaves nothing behind that changes them.

    Both run serially: with several OpenMP threads the powder sums come in a
    different order from run to run (last digits differ even between two
    fresh processes). 'set parallel' is a setting and outlasts reset, like
    variables do.

    The diffuse functions set every menu option they use, so this also
    passes without reset() (checked): nothing left over from the MMC run
    changes their results. reset() matters for raw run_macro() use."""
    import shutil
    import subprocess
    import conftest
    examples = REPO / "python" / "examples"
    script = tmp_path / "fresh.py"
    script.write_text(RESET_SCRIPT)
    env = os.environ | {"MPLBACKEND": "Agg"}
    if conftest.USE_SOURCE:
        env["PYTHONPATH"] = str(REPO / "python")
    r = subprocess.run([sys.executable, str(script), str(examples), str(tmp_path / "fresh.npz")],
                       env=env, capture_output=True, text=True, timeout=300)
    assert r.returncode == 0, r.stdout + r.stderr
    fresh = np.load(tmp_path / "fresh.npz")

    sro = tmp_path / "sro3"
    shutil.copytree(REPO / "tests" / "04.SRO" / "SRO_3", sro)
    text = (sro / "sro.3.mac").read_text()
    (sro / "sro_mc.mac").write_text(text[:text.index("@fourier.mac")] + "exit\n")
    dc.run_macro("set wait, off\nseed 1234\n@sro_mc.mac 0.5, 0.5, 0.5, 0.48, 0.48, 1.5",
                 cwd=sro)
    dc.reset()
    try:
        dc.run_macro("set parallel, useomp:serial")
        layer = dc.read_structure(examples / "cuau_layer.stru")
        f = dc.fourier(layer, ll=(0, 0, 0), lr=(2, 0, 0), ul=(0, 2, 0), npoints=(41, 41))
        ni = dc.read_structure(examples / "nickel.cell", supercell=3)
        p = dc.powder(ni, qmax=6.0)
    finally:
        dc.run_macro("set parallel, useomp:use, nthreads:all")     # the start-up default
    np.testing.assert_array_equal(f.intensity, fresh["fourier"])
    np.testing.assert_array_equal(p.intensity, fresh["powder"])


# ------------------------------------------------------- the C library
def test_backend(dc):
    assert dc.session().backend == "ctypes"


def test_value_precision(dc):
    """dc_value returns doubles: about 12 significant digits, what the
    expression evaluator keeps (the f2py module gave single precision)."""
    assert abs(dc.session().value("1.0/3.0") - 1.0 / 3.0) < 1e-11


def test_missing_library(monkeypatch):
    """Without the library, load() says where it looked."""
    from diffuse import _core
    monkeypatch.setattr(_core, "find_library", lambda: None)
    with pytest.raises(ImportError, match="libdiffuse_core"):
        _core.load()
