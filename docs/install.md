# Installation

There are no ready-made packages yet (wheels and a conda-forge package come
later in WP1). You build the package from the DiffuseCode sources with
`pip`, which compiles DISCUS on the way. That needs a Fortran compiler, CMake,
FFTW, HDF5 (with its Fortran interface) and readline. Two ways to get them:

* **System packages** (needs admin rights): see below.
* **Conda** (no admin rights): see [Without admin rights](#without-admin-rights-conda).

Graphics libraries (PGPLOT, X11) are not needed: the Python package is built
without them.

## With system packages (Debian, Ubuntu, WSL2)

Install the build tools and libraries:

```bash
sudo apt-get install -y gfortran cmake make libfftw3-dev libhdf5-dev \
    libreadline-dev zlib1g-dev python3-dev python3-venv
```

Then create a virtual environment called `.diffuse` in your home folder,
activate it, and install the package into it from the DiffuseCode folder:

```bash
python3 -m venv ~/.diffuse
source ~/.diffuse/bin/activate
pip install --upgrade pip
pip install .
```

Run `source ~/.diffuse/bin/activate` again in every new terminal before using
the package (or add it to `~/.bashrc`). Keeping the environment outside the
DiffuseCode folder means a fresh clone or `git clean` does not touch it.

This also installs NumPy, matplotlib and py3Dmol (3D views of structures). Use `.[test]` for what the tests
need, `.[docs]` for the documentation tools. The build takes a few minutes. Set
`CMAKE_BUILD_PARALLEL_LEVEL` to the number of CPU cores to make it faster,
for example `CMAKE_BUILD_PARALLEL_LEVEL=4 pip install .`.

Check that it works:

```bash
python python/examples/powder_demo.py
```

```{note}
Always install into a virtual environment (or a conda environment). A NumPy
installed with `pip install --user` next to the system's NumPy is what broke
the f2py step on Ubuntu 22.04 before.
```

## WSL2

WSL2 (Windows Subsystem for Linux) runs a real Linux, so the instructions
above work unchanged in its Ubuntu. Ubuntu 22.04 with Python 3.10 is tested
in CI. Three tips:

* **Work in the Linux file system** (`~/...`), not under `/mnt/c/...`: files
  on the Windows side are much slower to read and write from Linux, and the
  build reads and writes many files.
* **No display is needed.** The package has no graphics; matplotlib plots go
  to files, or to a window if your WSL2 has WSLg (Windows 11).
* **Jupyter:** see [Jupyter](#jupyter) below, and the next section for
  opening it in your Windows browser.

### JupyterLab on WSL2

JupyterLab runs inside WSL2; you use it from your Windows browser. Start it
without a browser on the Linux side:

```bash
jupyter lab --no-browser
```

It prints a link like `http://localhost:8888/lab?token=...`. Paste it into
the Windows browser (in Windows Terminal, Ctrl+click opens it): WSL2
forwards `localhost` to Windows.

To make that the default, so that a plain `jupyter lab` never tries to open
a browser, create the configuration file once:

```bash
jupyter lab --generate-config
```

and add this line to `~/.jupyter/jupyter_lab_config.py`:

```python
c.ServerApp.open_browser = False
```

If `localhost` does not reach JupyterLab from Windows (some WSL network
settings), let it listen on all interfaces and use the WSL address instead
of `localhost` in the link (the first address that `hostname -I` prints):

```bash
jupyter lab --no-browser --ip=0.0.0.0
```

Alternatively, let `jupyter lab` open the Windows browser itself: install
`wslu` (`sudo apt install wslu`), which provides `wslview`, and put these two
lines into `~/.jupyter/jupyter_lab_config.py` instead:

```python
c.ServerApp.browser = "wslview %s"
c.ServerApp.use_redirect_file = False
```

The second line is needed because by default JupyterLab opens a start-up
file whose path exists only inside WSL, which the Windows browser cannot
open; with it, the browser gets the link with the token.

## Jupyter

To use the package in JupyterLab with interactive plots (zoom, pan, a
readout of the coordinates under the mouse), install it with the `jupyter`
extra, which adds JupyterLab and [ipympl](https://matplotlib.org/ipympl/):

```bash
pip install ".[jupyter]"
jupyter lab
```

(The quotes stop the shell from reading the square brackets.) In the
notebook, before plotting:

```python
%matplotlib widget
```

Each figure then has a toolbar: the magnifier zooms into a box, the arrows
pan, home returns to the full view, the disk saves a PNG. Restart JupyterLab
after installing ipympl, or the figures show `Error displaying widget`.

If the notebook reports `No module named 'ipympl'` (or `'diffuse'`) although
you installed it, the notebook runs another Python than the environment you
installed into, usually because JupyterLab was started from elsewhere (an
older JupyterLab, or one that is still running). Check in a cell:

```python
import sys; print(sys.executable)
```

It should be the Python of your environment (e.g. `~/.diffuse/bin/python`).
Either start JupyterLab from that environment (`which jupyter` should point
into it), or install into the Python the notebook runs with `%pip install
ipympl` in a cell; then restart the kernel and reload the page.

Structures can be viewed in 3D in the notebook with
{meth}`Crystal.view() <diffuse.Crystal.view>` (py3Dmol, installed with the
package; no extension needed): rotate with the mouse, zoom with the wheel.

The notebook `python/examples/interactive_plots.ipynb` shows both with a
structure view, a powder pattern, a diffuse-scattering map and a PDF. Run a
copy of the examples folder, not the notebook in your DiffuseCode clone:
Jupyter saves a notebook with its results, so the file in the clone would
change and `git pull` would refuse to update it.

```bash
cp -r python/examples ~/diffuse-examples
cd ~/diffuse-examples
jupyter lab
```

(If it has happened: `git checkout -- python/examples/interactive_plots.ipynb`
throws away the saved results, then `git pull` works.)

## Without admin rights (conda)

Everything, including the compilers, can come from
[conda-forge](https://conda-forge.org). Install
[Miniforge](https://github.com/conda-forge/miniforge) into your home folder
(no admin rights needed), then, in the DiffuseCode folder:

```bash
conda env create -f environment.yml
conda activate diffuse
pip install --no-build-isolation .
```

`--no-build-isolation` makes pip build with the compilers and libraries of
the conda environment. This works on Linux and on HPC clusters; it is tested
in CI.

On a cluster you can also use its environment modules instead (for example
`module load gcc cmake hdf5 fftw`), then follow the system-package route from
`python3 -m venv` on.

## For development: a CMake build folder

To work on DiffuseCode itself, build with CMake directly and point Python at
the build folder and the source:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DDIFFUSE_PYTHON=ON -DDIFFUSE_GRAPHICS=OFF
cmake --build build -j
export DIFFUSE_PYTHON_PATH=$PWD/build/python
export PYTHONPATH=$PWD/python
```

The classic program without Python is described in {doc}`standalone`.

`DIFFUSE_PYTHON_PATH` takes precedence over an installed package, so the
tests (`pytest tests/python`) and the documentation then use the build folder
and the source in `python/`. Leave out `-DDIFFUSE_GRAPHICS=OFF`, and install
`libpng-dev libx11-dev pgplot5`, for the classic `discus_suite` with KUPLOT
plots.

## How it works, and its limits

This version drives DISCUS through the f2py module `suite_python`: Python
writes the commands to a small macro, DISCUS runs it in the same process, and
the results come back through files or KUPLOT data sets.

* DISCUS holds **one structure at a time** per process. A {class}`~diffuse.Crystal`
  remembers how it was made, and {func}`~diffuse.powder` reads it in again if
  another structure was loaded since.
* DISCUS output is captured, not printed. Set `DIFFUSE_ECHO=1`, or call
  `dc.session(echo=True)`, to see it as it is written, also in a Jupyter
  cell (useful for long calculations); the last output is in
  `dc.session().last_output`.
* A failed DISCUS command raises {class}`~diffuse.MacroError` with DISCUS's
  own error lines.
* {func}`~diffuse.run_macro` runs any DISCUS commands that are not wrapped yet.
