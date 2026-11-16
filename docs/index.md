# DiffuseCode from Python

```{warning}
This is an early preview. The Python interface covers a small set of DISCUS
commands so far, and the API may still change. It is built from the
`claude-python-app` branch.
```

The `diffuse` package runs [DISCUS](https://github.com/tproffen/DiffuseCode)
from Python. You read a structure, run a calculation, and get the results back
as NumPy arrays, ready for matplotlib, SciPy or your own analysis.

```python
import diffuse as dc

cr = dc.read_structure("anatase.cif", supercell=(10, 10, 4))
p = dc.powder(cr, qmin=0.5, qmax=10, dq=0.005, radiation="xray")
p.plot()
```

```{toctree}
:maxdepth: 2
:caption: Using it

install
gallery/index
api
```

```{toctree}
:maxdepth: 1
:caption: Development

_generated/testing
_generated/migration/DiffuseCode_Python_Migration_Plan
_generated/migration/WP0_Test_Baseline
_generated/migration/WP0_Test_Coverage_Gaps
```
