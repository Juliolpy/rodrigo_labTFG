# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
import os
import sys
sys.path.insert(0, os.path.abspath('..'))  # Añade el directorio del proyecto al path

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Rodrigo_lab'
copyright = '2025, Julio'
author = 'Julio'
release = '0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',      # Extrae docstrings de tu código Python
    'sphinx.ext.napoleon',     # Soporta Google style y NumPy style docstrings
    'sphinx.ext.viewcode',     # Muestra enlaces al código fuente
    'sphinx.ext.todo',         # Permite usar .. todo::
    'sphinx.ext.autosummary',  # Genera resúmenes automáticos de módulos y funciones
]
# Habilita la generación de archivos .rst automáticamente
autosummary_generate = True
# -- Configuración para los TODOs -------------------------------------------------

todo_include_todos = True

templates_path = ['_templates']
exclude_patterns = []

language = 'es'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']

# -- Autodoc: orden de miembros y estilo ------------------------------------------

autodoc_member_order = 'bysource'
autoclass_content = 'both'  # Incluye __init__ y clase docstring
