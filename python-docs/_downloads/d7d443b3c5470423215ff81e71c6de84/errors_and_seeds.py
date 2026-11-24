"""
Errors and reproducible random numbers
======================================

A failed DISCUS command raises :class:`diffuse.MacroError`, which carries
DISCUS's error number, type and message. The Python process keeps running
and the next call works.

:func:`diffuse.seed` makes random numbers repeat, and
:meth:`diffuse.Session.random_state` saves the generator's state so that a
calculation can be run again from the same point.
"""

import diffuse as dc

# %%
# Errors become exceptions
# ------------------------
try:
    dc.run_macro("discus\n  read\n    stru no_such_file.stru\n  exit\nexit")
except dc.MacroError as err:
    print("code:    ", err.code)
    print("category:", err.category)
    print("text:    ", err.text)
    print("details: ", err.details)

# %%
# The session is still usable after the error:
print("1 + 1 =", dc.session().value("1 + 1"))

# %%
# ``backend`` says how Python reaches DISCUS: ``"ctypes"`` (the C library
# libdiffuse_core, the default) or ``"f2py"`` (the older module).
print("backend:", dc.session().backend)

# %%
# Reproducible random numbers
# ---------------------------
# ``ran(0)`` is DISCUS's uniform random number between 0 and 1.
s = dc.session()

dc.seed(1234)
first = [s.value("ran(0)") for _ in range(3)]
dc.seed(1234)
again = [s.value("ran(0)") for _ in range(3)]
print(first)
print("same after seeding again:", first == again)

# %%
# Save the generator state, draw numbers, then go back:
state = s.random_state()
before = [s.value("ran(0)") for _ in range(3)]
s.set_random_state(state)
after = [s.value("ran(0)") for _ in range(3)]
print("same after restoring the state:", before == after)
