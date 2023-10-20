# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'MPCD'
copyright = '2023, Shendruk Lab'
author = 'Shendruk Lab'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx_toolbox.collapse']

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output
# "furo","alabastar","sphinx_rtd_theme"
html_theme = 'furo'
html_static_path = ['_static']
# html_theme_options = {
#     "light_css_variables": {
#         "color-brand-primary": "red",
#         "color-brand-content": "#CC3333",
#         "color-admonition-background": "orange",
#     },
# }
