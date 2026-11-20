# The standalone program (discus_suite)

The classic command-line program `discus_suite` (DISCUS, KUPLOT with
plots, DIFFEV, REFINE) is built from the same sources as the Python
package, and the migration keeps it working. CI builds it with graphics and
runs the legacy macro tests through it whenever the Fortran code changes
(see {doc}`_generated/testing`).

The Python work does not change how the program behaves. The library mode
used by Python (an error never ends the program) is switched on only by the
Python module.

## Dependencies (Debian, Ubuntu, WSL2)

```bash
sudo apt-get install -y gfortran cmake make libfftw3-dev libhdf5-dev \
    libreadline-dev zlib1g-dev libpng-dev libx11-dev pgplot5
```

`libpng-dev libx11-dev pgplot5` are for KUPLOT's plots; leave them out for
a build without graphics (below). CMake finds the Fortran compiler and
PGPLOT itself: no `FC` or `PGPLOT_DIR` is needed.

## Build

In the DiffuseCode folder:

```bash
cmake -S . -B build-standalone -DCMAKE_BUILD_TYPE=Release
cmake --build build-standalone -j
```

The program is `build-standalone/suite/prog/discus_suite`. It runs from
there, but it looks for its help files and manuals in a `share/` folder
next to its `bin/` folder (falling back to `/usr/local/share/`), so for
`help` and `manual` install it.

Use a build folder of its own: `pip install .` uses `build/`.

## Install

Into `/usr/local` (needs admin rights):

```bash
sudo cmake --install build-standalone
```

Without admin rights, choose a folder in your home when configuring, and
install there:

```bash
cmake -S . -B build-standalone -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX=$HOME/discus-suite
cmake --build build-standalone -j
cmake --install build-standalone
```

Then add `$HOME/discus-suite/bin` to your `PATH` (for example in
`~/.bashrc`). The install holds `bin/discus_suite` (and `pgxwin_server` for
the plots) and `share/` with the help files and the PDF manuals.

## Options

| CMake option | Default | Effect |
|---|---|---|
| `-DDIFFUSE_GRAPHICS=OFF` | `ON` | no PGPLOT, X11 or PNG: no plots on screen; for clusters and headless machines |
| `-DDISCUS_OMP=OFF` | `ON` | without OpenMP |
| `-DDIFFUSE_PYTHON=ON` | `OFF` | also build the Python module next to the program (development; see {doc}`install`) |
| `-DDEBUG=ON` | `OFF` | debug build with run-time checks |

## Updating

After `git pull`, build again; only the changed files are compiled:

```bash
cmake --build build-standalone -j
cmake --install build-standalone      # if you installed it
```

## Checking the build

The legacy test runner runs the test macros in `tests/` through the
program (a few minutes):

```bash
python3 tools/run_legacy_tests.py --exe build-standalone/suite/prog/discus_suite -j 4 --tags fast,mc
```

The installation guides in `Manual/` (`AAA_INSTALL_DISCUS_LINUX.pdf` and
others) describe the older procedure; where they differ, follow this page.
CI checks it on every Fortran change: it builds the program as above (with
the Python module added), installs it with `cmake --install`, and runs the
installed `discus_suite` on a small macro.
