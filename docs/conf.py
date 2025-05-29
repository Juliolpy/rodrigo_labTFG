import os
import sys
# 1) Ajusta este path al directorio que contenga tu paquete columbo_design
sys.path.insert(0, os.path.abspath('../src'))

# -- Project information -----------------------------------------------------
project = 'Columbo Design'
copyright = '2025, Julio'
author = 'Julio'

# -- General configuration ---------------------------------------------------
# Añade extensions para autodoc y, opcionalmente, napoleon si usas docstrings Google/NumPy
extensions = [
    'sphinx.ext.autodoc',
    # 'sphinx.ext.napoleon',  # descomenta si prefieres docstrings Google o NumPy
]

templates_path = ['_templates']
exclude_patterns = []

language = 'en'

# -- Options for HTML output -------------------------------------------------
html_theme = 'alabaster'
html_static_path = ['_static']
