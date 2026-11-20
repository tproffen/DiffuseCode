# DiffuseCode from Python

`python/diffuse/` is an early Python API for DISCUS. It currently covers a
small set of commands: reading a structure, and calculating powder patterns,
single-crystal intensities (including diffuse scattering) and pair
distribution functions.
The results are NumPy arrays, so they can be plotted with matplotlib or any
other library.

```python
import diffuse as dc

cr = dc.read_structure("examples/anatase.cif", supercell=(10, 10, 4))
print(cr)            # Crystal('anatase.cif', supercell=(10, 10, 4), natoms=4800, ...)

p = dc.powder(cr, qmin=0.5, qmax=10, dq=0.005, radiation="xray")
p.q, p.intensity     # NumPy arrays
p.plot()             # matplotlib; p.plot(x="tth") for 2-theta
```

The documentation (installation, API reference and an example gallery with
plots) is built from `docs/`; see "Documentation" below.

## Running the demo

Build DiffuseCode with the Python module (from the repository root). The
Python interface needs no graphics, so this builds without PGPLOT and X11
(see `docs/install.md` for the packages):

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DDIFFUSE_PYTHON=ON -DDIFFUSE_GRAPHICS=OFF
```

```bash
cmake --build build -j
```

Then run the demo (needs NumPy and matplotlib):

```bash
DIFFUSE_PYTHON_PATH=build/python PYTHONPATH=python python3 python/examples/powder_demo.py
```

It reads `examples/anatase.cif`, calculates X-ray and neutron patterns for a
10×10×4 supercell and an X-ray pattern for a 3×3×1 crystal, and shows the
plot.

For Jupyter, set the two variables before starting it, or at the top of the
notebook before `import diffuse`:

```python
import os, sys
os.environ["DIFFUSE_PYTHON_PATH"] = "/path/to/build/python"
sys.path.insert(0, "/path/to/DiffuseCode/python")
import diffuse as dc
```

## What is in the API

| Call | What it does |
|---|---|
| `read_structure(path, supercell=None)` | Reads `.cif` or `.cell` (asymmetric unit; space group applied, `supercell` cells built) or `.stru` (whole crystal, read as is). Returns a `Crystal`. |
| `Crystal.natoms`, `.lattice`, `.volume` | Atom count, `(a, b, c, alpha, beta, gamma)`, cell volume |
| `Crystal.atoms()`, `.positions()` | Atom names and fractional coordinates as arrays |
| `powder(crystal, qmin, qmax, dq, radiation, wavelength, method, temperature, lp_correction)` | Powder pattern on a Q grid. `radiation`: `"xray"`, `"neutron"`, `"electron"`; `method`: `"debye"` or `"complete"`. Returns a `PowderPattern`. |
| `PowderPattern.q`, `.intensity`, `.two_theta()`, `.plot()` | Results and a quick matplotlib plot |
| `fourier(crystal, ll, lr, ul=None, npoints, radiation, wavelength, temperature, subtract_average)` | Single-crystal intensities on a line (`ul` left out) or plane in reciprocal space; `subtract_average=100` leaves the diffuse scattering. Returns a `FourierMap`. |
| `FourierMap.intensity`, `.hkl()`, `.x`, `.y`, `.plot(log=...)` | Intensities (`[i, j]`: i along `ll`→`lr`, j along `ll`→`ul`), their (h, k, l), and an image or line plot |
| `pdf(crystal, rmax, dr, radiation, qmax, qdamp, qbroad, thermal, corrlin, corrquad, particle_diameter)` | Pair distribution function G(r), for a periodic crystal or a finite particle. Returns a `PDFPattern`. |
| `PDFPattern.r`, `.G`, `.plot()` | Results and a plot |
| `run_macro(text)` | Runs any DISCUS macro text, for commands that are not wrapped yet |
| `MacroError` | Raised when a DISCUS command fails; the message holds the `***` error lines |

## How it works, and its limits

This version drives the existing f2py module `suite_python`
(`python/suite_python.f90`). Python writes the commands to a small macro,
DISCUS runs it in the same process, and the results come back through
files or KUPLOT data sets. See `diffuse/_session.py`.

* DISCUS holds **one structure at a time** per process. A `Crystal` remembers
  how it was made, and `powder()` reads it in again if another structure was
  loaded since.
* DISCUS output is captured, not printed. Set `DIFFUSE_ECHO=1`, or call
  `dc.session(echo=True)`, to see it. The last output is in
  `dc.session().last_output`.
* Input files are copied to a scratch folder first, because DISCUS writes
  converted files (for example `.stru` from `.cif`) next to its input.

The migration plan (`migration/DiffuseCode_Python_Migration_Plan.md`,
WP1–WP5) replaces `suite_python` underneath with a `bind(C)` core and a
proper `pip install`. The aim is to keep this API while doing so.

## Documentation

Sphinx sources are in `docs/`: Markdown pages, the API reference generated
from the docstrings (NumPy style), and a gallery made by running every script
in `python/examples/`. To build it:

```bash
pip install -r docs/requirements.txt
```

```bash
DIFFUSE_PYTHON_PATH=build/python python3 -m sphinx -W -b html docs docs/_build/html
```

`DIFFUSE_DOCS_NO_GALLERY=1` skips the gallery (no compiled module needed);
`DIFFUSE_DOCS_OFFLINE=1` skips the links into the NumPy/matplotlib docs.
CI builds the docs on every push and keeps the HTML as the `docs-html`
artifact of the run. Pushes to `claude-python-app` also publish them at
https://tproffen.github.io/DiffuseCode/python-docs/. The `gh-pages` branch
is the DISCUS homepage, so `tools/publish_docs.sh` only replaces the
`python-docs/` folder (plus an empty `.nojekyll` at the root, which Sphinx
output needs) and refuses to publish if anything else would change. The
docs are not linked from the homepage.

When adding to the API: give every public function a NumPy-style docstring
(Parameters, Returns, Raises, Examples), list it in `docs/api.md`, and add an
example script to `python/examples/` for each new capability (it becomes a
gallery page and is run by the tests).

## Tests

`tests/python/test_diffuse_api.py` covers this API, and CI runs it together
with the demo script. See `tests/README.md`.
