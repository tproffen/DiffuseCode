# Installation

There is no `pip install` yet (that comes with the packaging work, WP1 in
the migration plan). For now, build DiffuseCode with its Python module and
point Python at the build.

## Build

On Debian or Ubuntu, install the build dependencies:

```bash
sudo apt-get install -y gfortran cmake libfftw3-dev libhdf5-dev libreadline-dev \
    libpng-dev libx11-dev zlib1g-dev pgplot5 python3-dev python3-numpy python3-matplotlib
```

Then, from the repository root:

```bash
FC=gfortran CC=gcc CXX=g++ cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
    -DDIFFUSE_PYTHON=ON -DPGPLOT_DIR=/usr/lib/pgplot5
cmake --build build -j
```

This builds `discus_suite` and the Python module `build/python/suite_python*.so`.

## Use

Two settings tell Python where things are:

* `DIFFUSE_PYTHON_PATH`: the folder with the compiled module (`build/python`).
* the `python/` folder of the repository on the Python path, for the
  `diffuse` package.

```bash
export DIFFUSE_PYTHON_PATH=$PWD/build/python
export PYTHONPATH=$PWD/python
python3 python/examples/powder_demo.py
```

In a Jupyter notebook, set them before `import diffuse`:

```python
import os, sys
os.environ["DIFFUSE_PYTHON_PATH"] = "/path/to/DiffuseCode/build/python"
sys.path.insert(0, "/path/to/DiffuseCode/python")
import diffuse as dc
```

## How it works, and its limits

This version drives DISCUS through the f2py module `suite_python`: Python
writes the commands to a small macro, DISCUS runs it in the same process, and
the results come back through files or KUPLOT data sets.

* DISCUS holds **one structure at a time** per process. A {class}`~diffuse.Crystal`
  remembers how it was made, and {func}`~diffuse.powder` reads it in again if
  another structure was loaded since.
* DISCUS output is captured, not printed. Set `DIFFUSE_ECHO=1`, or call
  `dc.session(echo=True)`, to see it; the last output is in
  `dc.session().last_output`.
* A failed DISCUS command raises {class}`~diffuse.MacroError` with DISCUS's
  own error lines.
* {func}`~diffuse.run_macro` runs any DISCUS commands that are not wrapped yet.
