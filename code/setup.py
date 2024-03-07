from setuptools import Extension, setup
from Cython.Build import cythonize

compiler_directives = {}
compiler_directives['language_level'] = 3
define_macros = []

debug = False
if debug:
    define_macros.append(('CYTHON_TRACE', '1'))
    compiler_directives['profile'] = True
    compiler_directives['linetrace'] = True

ext_modules=[
    Extension("grid", ["grid.pyx"], include_dirs = ['.'], define_macros = define_macros),
    Extension("cellDynamics", ["cellDynamics.pyx"], include_dirs = ['.'], define_macros = define_macros),
    Extension("spatialDynamics", ["spatialDynamics.pyx"], include_dirs = ['.'], define_macros = define_macros),
    Extension("wind", ["wind.pyx"], include_dirs = ['.'], define_macros = define_macros),
    Extension("spatialDependence", ["spatialDependence.pyx"], include_dirs = ['.'], define_macros = define_macros),
    Extension("populationsAndParameters", ["populationsAndParameters.pyx"], include_dirs = ['.'], define_macros = define_macros)
]

setup(ext_modules = cythonize(
    ext_modules,
    compiler_directives = compiler_directives,
    annotate = True,
    force = True))

