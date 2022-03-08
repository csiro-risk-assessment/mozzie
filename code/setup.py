import os
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules=[
    Extension("grid", ["grid.pyx"], include_dirs=['.']),
    Extension("cellDynamics", ["cellDynamics.pyx"], include_dirs=['.']),
    Extension("spatialDynamics", ["spatialDynamics.pyx"], include_dirs=['.']),
    Extension("wind", ["wind.pyx"], include_dirs=['.']),
    Extension("spatialDependence", ["spatialDependence.pyx"], include_dirs=['.']),
    Extension("populationsAndParameters", ["populationsAndParameters.pyx"], include_dirs=['.'])
]

setup(
  name="mozzie",
  ext_modules=ext_modules,
  cmdclass = {'build_ext': build_ext},
  script_args = ['build_ext'],
  options = {'build_ext':{'inplace':True, 'force':True}}
 )

