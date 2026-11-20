"""WP2 core: library mode, errors as exceptions, reset and random state.

The "done when" test of WP2 in the migration plan: a deliberately bad input
raises a Python exception, the process survives, nothing reaches the
terminal, and the next call works.

    DIFFUSE_PYTHON_PATH=<build>/python pytest tests/python/test_wp2_core.py
"""

from __future__ import annotations

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
