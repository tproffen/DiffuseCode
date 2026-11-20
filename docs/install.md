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

This also installs NumPy and matplotlib. Use `.[test]` for what the tests
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
* **Jupyter:** run `pip install jupyterlab` in the virtual environment and
  start `jupyter lab`; open the link it prints in your Windows browser.

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
  `dc.session(echo=True)`, to see it; the last output is in
  `dc.session().last_output`.
* A failed DISCUS command raises {class}`~diffuse.MacroError` with DISCUS's
  own error lines.
* {func}`~diffuse.run_macro` runs any DISCUS commands that are not wrapped yet.
