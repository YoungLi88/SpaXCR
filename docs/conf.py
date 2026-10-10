from __future__ import annotations

import os
import sys
from pathlib import Path

DOCS_DIR = Path(__file__).resolve().parent
ROOT_DIR = DOCS_DIR.parent

sys.path.insert(0, str(ROOT_DIR))

project = "spaxcr"
author = "spaxcr contributors"
copyright = "2026, spaxcr contributors"
version = "0.1.1"
release = "0.1.1"

extensions = [
    "myst_parser",
    "sphinx_copybutton",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

master_doc = "index"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "sphinx_rtd_theme"
html_title = "spaxcr documentation"
html_short_title = "spaxcr"
html_static_path = ["_static"]
html_css_files = ["custom.css"]

html_theme_options = {
    "collapse_navigation": False,
    "navigation_depth": 3,
    "sticky_navigation": True,
    "includehidden": True,
    "titles_only": False,
}

myst_heading_anchors = 3
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "substitution",
]

copybutton_prompt_text = r">>> |\.\.\. |\$ "
copybutton_prompt_is_regexp = True

nitpicky = False
suppress_warnings = [
    "myst.header",
]
