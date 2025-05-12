# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'scanpy_plus'
copyright = '2025, Vijaya Baskar, Masami Kuri'
author = 'Vijaya Baskar, Masami Kuri'
release = 'v1.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = []

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_theme = 'furo'
html_static_path = ['_static']

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",          # For Google/NumPy style docstrings
    "sphinx_autodoc_typehints",     # Optional: for type hints
]

# Optional: set the module import path
import os
import sys
sys.path.insert(0, os.path.abspath('../../'))  # Adjust if needed

#autodoc_default_options = {
#    'members': True,
#    'undoc-members': True,
#    'show-inheritance': True,
#}