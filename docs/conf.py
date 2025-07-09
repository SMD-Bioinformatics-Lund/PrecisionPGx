import os
import sys

# Add project root to sys.path to enable Sphinx to find modules
sys.path.insert(0, os.path.abspath('..'))

project = 'PrecisionPGx'
author = 'PrecisionPGx Team'

# The full version, including alpha/beta/rc tags
release = '0.1.0'

extensions = ['sphinx.ext.autodoc', 'sphinx.ext.napoleon']

templates_path = ['_templates']
exclude_patterns = []

html_theme = 'alabaster'
html_static_path = ['_static']