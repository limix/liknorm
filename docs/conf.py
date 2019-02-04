from __future__ import unicode_literals

import os
import subprocess
import sys

folder = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(folder, "..", "VERSION")) as f:
    version = f.read().strip()

extensions = [
    "sphinx.ext.viewcode",
    "breathe",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
]

breathe_projects = {"liknorm": "doxyxml/"}
breathe_default_project = "liknorm"
napoleon_numpy_docstring = True

project = "liknorm"
copyright = "2019, Danilo Horta"
author = "Danilo Horta"

release = version
templates_path = ["_templates"]
source_suffix = ".rst"
master_doc = "index"
language = "en"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "conf.py"]
pygments_style = "default"
todo_include_todos = False
primary_domain = "c"
highlight_language = "none"

html_static_path = ["_static"]
html_theme = "sphinx_rtd_theme"
html_sidebars = {
    "**": ["about.html", "searchbox.html", "localtoc.html", "relations.html"]
}
htmlhelp_basename = "liknormdoc"


def run_doxygen(folder):
    """Run the doxygen make command in the designated folder"""

    try:
        retcode = subprocess.call("cd %s; make doxygen" % folder, shell=True)
        if retcode < 0:
            sys.stderr.write("doxygen terminated by signal %s" % (-retcode))
    except OSError as e:
        sys.stderr.write("doxygen execution failed: %s" % e)


def generate_doxygen_xml(app):
    """Run the doxygen make commands if we're on the ReadTheDocs server"""

    read_the_docs_build = os.environ.get("READTHEDOCS", None) == "True"

    if read_the_docs_build:
        run_doxygen(".")


def setup(app):

    # Add hook for building doxygen xml when needed
    app.connect("builder-inited", generate_doxygen_xml)
