import pathlib
import os
import sys

# project root directory to the sys.path
project_root = pathlib.Path(__file__).parents[2].resolve().as_posix()
module_dir   = os.path.join(project_root,'geoteqpy')

sys.path.insert(0, project_root)
sys.path.insert(0, module_dir)

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'GeoTEQpy'
copyright = '2024, Anthony Jourdon'
author = 'Anthony Jourdon'
release = '1.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
  'sphinx.ext.duration',
  'sphinx.ext.intersphinx',
  'sphinx.ext.doctest',
  'sphinx.ext.autodoc',
  ]

templates_path = ['_templates']
exclude_patterns = []
include_patterns = ['**',
                    '../geoteqpy/**']


rst_prolog = """
.. _pTatin3d: https://github.com/laetitialp/ptatin-gene
.. _ASAGI: https://github.com/TUM-I5/ASAGI
.. _SeisSol: https://seissol.org/
"""

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

html_theme_options = {'prev_next_buttons_location': 'both',
                      'style_external_links': True,}
