from setuptools import Extension, setup
from Cython.Build import cythonize

compiler_directives = {}
compiler_directives['language_level'] = 3

compiler_directives['profile'] = True
compiler_directives['linetrace'] = True

define_macros = []
define_macros.append(('CYTHON_TRACE', '1'))

setup(ext_modules = cythonize(Extension(
       "cyApple",
       sources=["*.pyx"],
       define_macros=define_macros,
  ), compiler_directives=compiler_directives))
