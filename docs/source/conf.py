#====================================================================================================
#  Copyright (c) 2024, 
#  Anthony Jourdon, 
#
#  project:  GeoTEQpy
#  filename: conf.py
#
#  This file is part of GeoTEQpy.
#
#  GeoTEQpy is free software: you can redistribute it and/or modify it under the terms 
#  of the GNU General Public License as published by the Free Software Foundation, either 
#  version 3 of the License, or any later version.
#
#  GeoTEQpy is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
#  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#  See the GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License along with GeoTEQpy. 
#  If not, see <https://www.gnu.org/licenses/>.
#====================================================================================================

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
.. |paper-title| replace:: *"3D reconstruction of complex fault systems from volumetric geodynamic shear zones using medial axis transform"*
.. _repository: https://
"""

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = []

html_theme_options = {'prev_next_buttons_location': 'both',
                      'style_external_links': True,}
