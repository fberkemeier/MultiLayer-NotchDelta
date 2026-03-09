from datetime import datetime
import os
import sys

sys.path.insert(0, os.path.abspath("../../src"))

project = "MultiLayer-NotchDelta"
author = "Francisco Berkemeier"
copyright = f"{datetime.now().year}, {author}"
release = "0.1.0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_copybutton",
    "myst_parser",
]

templates_path = ["_templates"]
exclude_patterns = [".ipynb_checkpoints", "**/.ipynb_checkpoints"]

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

html_theme_options = {
    "version_selector": False,
    "language_selector": False,
    "flyout_display": "hidden",
}

