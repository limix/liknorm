import os
import sphinx_rtd_theme

folder = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(folder, '../VERSION')) as f:
    version = f.read().strip()

extensions = [
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode'
]
templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'
project = 'liknorm'
copyright = '2016, Danilo Horta'
author = 'Danilo Horta'
release = version
language = None
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
pygments_style = 'sphinx'
todo_include_todos = False
html_theme = 'sphinx_rtd_theme'
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
html_static_path = ['_static']
htmlhelp_basename = 'liknormdoc'
latex_elements = {}
latex_documents = [
    (master_doc, 'liknorm.tex', 'liknorm Documentation',
     'Danilo Horta', 'manual'),
]
man_pages = [
    (master_doc, 'liknorm', 'liknorm Documentation',
     [author], 1)
]
texinfo_documents = [
    (master_doc, 'liknorm', 'liknorm Documentation',
     author, 'liknorm', 'One line description of project.',
     'Miscellaneous'),
]
