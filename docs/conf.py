"""Sphinx configuration for the DiffuseCode Python documentation.

Build (from the repository root), with the package installed:

    pip install .[docs]
    python3 -m sphinx -W -b html docs docs/_build/html

or against a CMake build folder:

    pip install -r docs/requirements.txt
    DIFFUSE_PYTHON_PATH=build/python python3 -m sphinx -W -b html docs docs/_build/html

The example gallery runs python/examples/*.py and needs the compiled core
libdiffuse_core (DIFFUSE_PYTHON_PATH, or an installed package). Without it, set DIFFUSE_DOCS_NO_GALLERY=1 to
build everything except the gallery. Without network access, set
DIFFUSE_DOCS_OFFLINE=1 (no links into the Python/NumPy/matplotlib docs).
"""

import os
import shutil
import sys
from pathlib import Path

DOCS = Path(__file__).resolve().parent
REPO = DOCS.parent
# Document (and run the gallery with) the source in python/ when building
# against a CMake build folder (DIFFUSE_PYTHON_PATH) or when diffuse is not
# installed; otherwise the installed package ("pip install .").
import importlib.util  # noqa: E402

if os.environ.get("DIFFUSE_PYTHON_PATH") or importlib.util.find_spec("diffuse") is None:
    sys.path.insert(0, str(REPO / "python"))

project = "DiffuseCode Python"
author = "The DiffuseCode developers"
copyright = "DiffuseCode developers"

extensions = [
    "myst_parser",               # Markdown pages
    "sphinx.ext.autodoc",        # API pages from docstrings
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",       # NumPy-style docstrings
    "sphinx.ext.intersphinx",    # links to the NumPy, matplotlib, Python docs
    "sphinx.ext.viewcode",
]

source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
exclude_patterns = ["_build", "**.ipynb_checkpoints"]
myst_heading_anchors = 3

# -- API docs -----------------------------------------------------------------
autodoc_member_order = "bysource"
autodoc_typehints = "description"
autodoc_default_options = {"members": True}
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_rtype = False
# diffuse loads libdiffuse_core only when a session starts, so importing it
# for the API pages works without the compiled library.

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
}
if os.environ.get("DIFFUSE_DOCS_OFFLINE"):     # no network: no outside links
    intersphinx_mapping = {}

# The migration documents skip heading levels in places; that is fine there.
suppress_warnings = ["myst.header"]

# -- Pages taken from elsewhere in the repository -------------------------------
# The migration documents stay in migration/ (one source of truth); copy them
# into the docs tree so that the links between them keep working.
_generated = DOCS / "_generated"
shutil.rmtree(_generated, ignore_errors=True)
(_generated / "migration").mkdir(parents=True)
for md in sorted((REPO / "migration").glob("*.md")):
    shutil.copy(md, _generated / "migration" / md.name)
shutil.copy(REPO / "tests" / "README.md", _generated / "testing.md")

# -- Example gallery ------------------------------------------------------------
if not os.environ.get("DIFFUSE_DOCS_NO_GALLERY"):
    extensions.append("sphinx_gallery.gen_gallery")
    sphinx_gallery_conf = {
        "examples_dirs": str(REPO / "python" / "examples"),
        "gallery_dirs": "gallery",
        "filename_pattern": r"\.py$",      # run every example
        "download_all_examples": False,
        "matplotlib_animations": False,
        "image_scrapers": ("matplotlib",),
    }
else:
    # Keep the table of contents valid: a placeholder instead of the gallery.
    shutil.rmtree(DOCS / "gallery", ignore_errors=True)
    (DOCS / "gallery").mkdir()
    (DOCS / "gallery" / "index.md").write_text(
        "# Examples\n\nThe example gallery was not built "
        "(`DIFFUSE_DOCS_NO_GALLERY` was set).\n")
    exclude_patterns.append("sg_execution_times.rst")

# -- HTML -------------------------------------------------------------------------
html_theme = "pydata_sphinx_theme"
html_title = "DiffuseCode Python"
html_theme_options = {
    "github_url": "https://github.com/tproffen/DiffuseCode",
    "navigation_with_keys": False,
    "show_toc_level": 2,
    "header_links_before_dropdown": 3,
}
html_static_path = []
