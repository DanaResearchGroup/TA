# Sphinx configuration for the TA documentation.

project = "TA"
copyright = "2026, Alon Grinberg Dana"
author = "Alon Grinberg Dana"
release = "0.1.0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx.ext.mathjax",
    "myst_parser",
]

# MyST (Markdown) settings
myst_enable_extensions = [
    "colon_fence",
    "dollarmath",
]
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# Theme
html_theme = "furo"
html_static_path = ["_static"]
html_title = "TA Documentation"

# Autodoc
autodoc_member_order = "bysource"
autodoc_typehints = "description"

# Napoleon (Google/NumPy style docstrings)
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_rtype = False

# Intersphinx
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "cantera": ("https://cantera.org/documentation/docs-3.0/sphinx/html/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
}
