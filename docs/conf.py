#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# License           : GPL3
# Author            : Jingxin Fu <jingxinfu.tj@gmail.com>
# Date              : 10/02/2020
# Last Modified Date: 11/02/2020
# Last Modified By  : Jingxin Fu <jingxinfu.tj@gmail.com>
# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('..'))


# -- Project information -----------------------------------------------------

project = 'Biopyutils'
import time
copyright = u'2020-{}'.format(time.strftime("%Y"))
author = 'Jingxin Fu'
import Biopyutils
release = Biopyutils.__version__

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.githubpages',
    'sphinx.ext.doctest',
    'sphinx.ext.intersphinx',
    'sphinx.ext.autosummary',
    'numpydoc',
]
autosummary_generate = True
numpydoc_show_class_members = False
# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

source_suffix = '.rst'
# The master toctree document.
master_doc = 'index'

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = 'bootstrap'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_theme_options = {
            'source_link_position': "footer",
            'bootswatch_theme': "paper",
            'navbar_sidebarrel': False,
            'bootstrap_version': "3",
        }
# Add any paths that contain custom themes here, relative to this directory.
import sphinx_bootstrap_theme
html_theme_path = sphinx_bootstrap_theme.get_html_theme_path()

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
html_favicon = "_static/favicon.ico"

# If true, links to the reST sources are added to the pages.
html_show_sourcelink = False
# Output file base name for HTML help builder.
htmlhelp_basename = 'Bipyutilsdoc'

# -- Extension configuration -------------------------------------------------
# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'
# -- Options for intersphinx extension ---------------------------------------

# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {'https://docs.python.org/3/': None}
