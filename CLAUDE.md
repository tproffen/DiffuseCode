# DiffuseCode – notes for Claude Code sessions

## Branches

Until the first release of the Python package, **all migration work is built
against and merged into `claude-python-app`**, not `master`. Start feature
branches from `origin/claude-python-app`.

When a piece of work is finished and pushed, **open the pull request
yourself, with base `claude-python-app`** (never `master`), so the owner
cannot merge into `master` by accident. Say in the PR description that it
targets `claude-python-app`.

## Plan and progress

- `migration/DiffuseCode_Python_Migration_Plan.md` – overall plan (option B,
  work packages WP0–WP10).
- `migration/WP0_Test_Baseline.md` – test baseline tasks and progress log.
- `tests/README.md` – how to run the tests.

## Build and test

Users: `pip install .` in a virtual environment (pyproject.toml,
scikit-build-core, graphics off); `environment.yml` for conda (no admin
rights). See `docs/install.md`. Development with a CMake build folder:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DDIFFUSE_PYTHON=ON
# add -DDIFFUSE_GRAPHICS=OFF for a build without PGPLOT/X11 (what the Python
# package will use; CI tests both)
cmake --build build -j
pip install -r tests/requirements.txt          # gemmi
DIFFUSE_PYTHON_PATH=build/python python3 -m pytest tests/python
python3 tools/run_legacy_tests.py --exe build/suite/prog/discus_suite -j 4 --tags fast,mc
```

CMake finds the compilers itself (`FC` is optional). Keep build trees outside `tests/`.

CI runs only the jobs a change needs (`tools/ci_changes.py`, see
`tests/README.md`). After substantial changes (Fortran core, build or
packaging, several areas at once), **offer the owner a full CI run** and
start it only after they say yes: `actions_run_trigger` with `run_workflow`,
`tests.yml`, the branch as ref. The "Run workflow" button does not exist
before the workflow is on `master` (first release).

## Documentation

Sphinx in `docs/` (MyST Markdown, autodoc + napoleon, sphinx-gallery over
`python/examples/`). Build with
`DIFFUSE_PYTHON_PATH=build/python python3 -m sphinx -W -b html docs docs/_build/html`
(`DIFFUSE_DOCS_OFFLINE=1` without network). New public API needs a NumPy-style
docstring, an entry in `docs/api.md`, and an example in `python/examples/`.
CI uploads the HTML as the `docs-html` artifact, and pushes to
`claude-python-app` publish it to
https://tproffen.github.io/DiffuseCode/python-docs/ (job `publish-docs`).
The `gh-pages` branch holds the DISCUS homepage: publishing may only change
`python-docs/` and the empty root `.nojekyll` (agreed with the owner, Nov
2026). `tools/publish_docs.sh` enforces this; keep it that way, and do not
link the docs from the homepage.
