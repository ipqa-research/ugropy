import os
import pathlib
import sys
import toml


CURRENT_PATH = pathlib.Path(os.path.abspath(os.path.dirname(__file__)))
UGROPY_PATH = CURRENT_PATH.parent.parent

sys.path.insert(0, str(UGROPY_PATH))

# Get release from pyproject
pyproject_path = UGROPY_PATH / "pyproject.toml"
with open(pyproject_path, "r") as f:
    pyproject_toml = toml.load(f)

project_version = pyproject_toml["project"]["version"]

project = "ugropy"
copyright = "2023, Salvador Eduardo Brandolín"
author = "Salvador Eduardo Brandolín"
release = project_version

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosummary",
    "nbsphinx",
    "sphinx_copybutton",
    "sphinxcontrib.bibtex",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "**.ipynb_checkpoints"]

# =============================================================================
# EXTRA CONF
# =============================================================================
# nbsphinx
nbsphinx_execute = "always"

autodoc_member_order = "bysource"

bibtex_bibfiles = ["refs.bib"]

add_module_names = False
# =============================================================================
# NUMPY DOC
# =============================================================================
numpydoc_class_members_toctree = False

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

# The master toctree document.
master_doc = "index"

"""Configuration Module."""
# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ["_static"]
# html_favicon = "_static/favicon.ico"
