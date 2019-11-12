#!/bin/bash
module load python/3.7.2

cythonize -3 -a *.pyx

# get the numpy include directory using import numpy then numpy.get_include()
gcc_flags="-shared -fno-strict-aliasing -Wsign-compare -Wunreachable-code -DNDEBUG -g -O3 -Wall -I/apps/python/3.7.2/include/python3.7m -I/apps/python/3.7.2/lib/python3.7/site-packages/numpy-1.15.4-py3.7-linux-x86_64.egg/numpy/core/include -fPIC"
mpicc $gcc_flags cellDynamics.c -o cellDynamics.so
mpicc $gcc_flags grid.c -o grid.so
mpicc $gcc_flags wind.c -o wind.so
mpicc $gcc_flags spatialDynamics.c -o spatialDynamics.so
mpicc $gcc_flags populationsAndParameters.c -o populationsAndParameters.so
mpicc $gcc_flags spatialDependence.c -o spatialDependence.so
