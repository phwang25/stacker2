# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sys
sys.path.insert(0, os.path.abspath('../'))

project = 'StACKER'
copyright = '2024, Eric D. Sakkas'
author = 'Eric D. Sakkas'
release = '1.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'numpydoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.todo',
    'sphinx.ext.viewcode',
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.coverage',
    'sphinx.ext.doctest',
    'sphinx.ext.graphviz',
    'sphinx.ext.ifconfig',
    'matplotlib.sphinxext.plot_directive',
    'sphinx_design',
    'sphinx.ext.imgmath',
    ]

autosummary_generate = True

napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_use_param = True
napoleon_use_ivar = True

toc_object_entries_show_parents = "hide"

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

autodoc_typehints = 'none'
autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'private-members': True,
    'special-members': '__init__',
    'inherited-members': True,
    'show-inheritance': True,
}

autodoc_member_order = 'bysource'
add_module_names = False

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output
pygments_style = 'sphinx'
html_theme = "sphinx_book_theme"
html_static_path = ['_static']
html_css_files = [
    'numpy.css',
]

html_logo = "images/stackerlogo.svg"

html_use_modindex = True
html_copy_source = False
html_domain_indices = False
html_file_suffix = '.html'
html_title = "StACKER Documentation"

html_theme_options = {
    "logo": {
        "image_light": "images/stackerlogo.svg",
        "image_dark": "images/stacker_logo_dark.png",
    },
    "home_page_in_toc": True,
    "icon_links" : [
        {
            "name" : "PyPI",
            "url" : "https://pypi.org/project/pistacker/",
            "type": "url",
            "icon" : "https://raw.githubusercontent.com/esakkas24/stacker/refs/heads/main/docs/images/pypi.svg"
        }
    ],
    "navbar_start" : [],
    "navbar_center": ["navbar-nav"],
    "navbar_end": [
        "navbar-icon-links",
        "search-button-field"
    ],
    "header_links_before_dropdown": 4,
    "navbar_persistent": [],
    "github_url" : "https://github.com/esakkas24/stacker",
    "repository_url" : "https://github.com/esakkas24/stacker",
    "use_repository_button" : True,
    "collapse_navigation": True,
}


