Diffuse program collection
================================================================

Welcome to the 'diffuse program' collection.
The website can be found at http://tproffen.github.io/DiffuseCode

**Documentation for the Python package and for building from source:**
https://tproffen.github.io/DiffuseCode/python-docs/

* [Installation of the Python package](https://tproffen.github.io/DiffuseCode/python-docs/install.html) (pip in a
  virtual environment, conda without admin rights, WSL2)
* [Building the standalone program discus_suite](https://tproffen.github.io/DiffuseCode/python-docs/standalone.html)
  (Linux, WSL2; install without admin rights)
* [Examples](https://tproffen.github.io/DiffuseCode/python-docs/gallery/index.html) and the
  [API reference](https://tproffen.github.io/DiffuseCode/python-docs/api.html)

Depending on the package you have downloaded you will find one or all of the
following program directories after unpacking the archive.

DISCUS_SUITE :  The Main program for DISCUS/KUPLOT/DIFFEV/REFINE
                This is the standard program that should be used
                and includes all other sections under a common roof top
AUTHOR  : R.B. Neder  (reinhard.neder@fau.de)

DISCUS  : Diffuse Scattering & Defect Structure Simulation
AUTHORS : R.B. Neder  (reinhard.neder@fau.de)
          Th. Proffen (tproffen@ornl.gov)

DIFFEV  : Generic refinement program using evolutionary algoritm
AUTHOR  : R.B. Neder  (reinhard.neder@fau.de)

KUPLOT  : General plotting program (well suited for DISCUS output)
AUTHOR  : Th. Proffen (tproffen@ornl.gov)
          R.B. Neder  (reinhard.neder@fau.de)

MIXSCAT : Program to generate differential PDFs from n/X data
AUTHOR  : C. Wurden K. Page A. Llobet
          Th. Proffen (tproffen@ornl.gov)

PYTHON PACKAGE (preview)
================================================================

The `diffuse` Python package runs DISCUS from Python and returns the results
as NumPy arrays. It is an early preview, developed on the
`claude-python-app` branch. After installing the system packages listed in
the [installation guide](https://tproffen.github.io/DiffuseCode/python-docs/install.html)
(which also has the conda route for users without admin rights), in the
DiffuseCode folder:

```bash
python3 -m venv ~/.diffuse
source ~/.diffuse/bin/activate
pip install .
```

INSTALLATION
================================================================

Binary
------

Binary distributions as well as the source code can be found 
at https://github.com/tproffen/DiffuseCode/releases

Linux (and Windows with WSL2)
-----

To build from source, follow [Building the standalone program](https://tproffen.github.io/DiffuseCode/python-docs/standalone.html)
(kept up to date and checked by CI). The older guide is
[AAA_INSTALL_DISCUS_LINUX.pdf](Manual/AAA_INSTALL_DISCUS_LINUX.pdf).

OSX - Using homebrew
---

Instructions to build from source please check file [AAA_INSTALL_DISCUS_MACOS.pdf](Manual/AAA_INSTALL_DISCUS_MACOS.pdf).

Windows
-------

Recommend you install the binary found [here](https://github.com/tproffen/DiffuseCode/releases).

CYGWIN
------

Instructions to build from source please check file [AAA_INSTALL_DISCUS_CYGWIN.pdf](Manual/AAA_INSTALL_DISCUS_CYGWIN.pdf).


[![tests](https://github.com/tproffen/DiffuseCode/actions/workflows/tests.yml/badge.svg?branch=claude-python-app)](https://github.com/tproffen/DiffuseCode/actions/workflows/tests.yml)
