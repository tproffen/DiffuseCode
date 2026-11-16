# DiffuseCode – notes for Claude Code sessions

## Branches

Until the first release of the Python package, **all migration work is built
against and merged into `claude-python-app`**, not `master`. Start feature
branches from `origin/claude-python-app` and target pull requests at it.

## Plan and progress

- `migration/DiffuseCode_Python_Migration_Plan.md` – overall plan (option B,
  work packages WP0–WP10).
- `migration/WP0_Test_Baseline.md` – test baseline tasks and progress log.
- `tests/README.md` – how to run the tests.

## Build and test

```bash
FC=gfortran CC=gcc CXX=g++ cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
    -DDIFFUSE_PYTHON=ON -DPGPLOT_DIR=/usr/lib/pgplot5
cmake --build build -j
DIFFUSE_PYTHON_PATH=build/python python3 -m pytest tests/python
python3 tools/run_legacy_tests.py --exe build/suite/prog/discus_suite -j 4 --tags fast,mc
```

CMake needs `FC` set. Keep build trees outside `tests/`.

## Documentation

Sphinx in `docs/` (MyST Markdown, autodoc + napoleon, sphinx-gallery over
`python/examples/`). Build with
`DIFFUSE_PYTHON_PATH=build/python python3 -m sphinx -W -b html docs docs/_build/html`
(`DIFFUSE_DOCS_OFFLINE=1` without network). New public API needs a NumPy-style
docstring, an entry in `docs/api.md`, and an example in `python/examples/`.
CI uploads the HTML as the `docs-html` artifact. Decision (Nov 2026): the
docs are **not published** until the first release; the artifact is the way
to read them. Do not publish to the `gh-pages` branch: it holds the DISCUS
homepage. When publishing is wanted, the preferred route is a separate docs
repository with its own GitHub Pages site.

