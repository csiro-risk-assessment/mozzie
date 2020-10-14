#!/bin/bash

module load openmpi/3.1.4-ofed45-gcc
module load gcc
module load git
module load python/3.7.2

cd auxillary; ./build_pearcey.sh ; cd ..

# get the numpy include directory using import numpy then numpy.get_include()
gcc_flags="-shared -fno-strict-aliasing -Wsign-compare -Wunreachable-code -DNDEBUG -g -O3 -Wall -I/apps/python/3.7.2/include/python3.7m -I/apps/python/3.7.2/lib/python3.7/site-packages/numpy-1.15.4-py3.7-linux-x86_64.egg/numpy/core/include -fPIC -march=native"

mpicc $gcc_flags csvparser.c -o csvparser.so

cythonize -3 -a *.pyx

mpicc $gcc_flags cellDynamics.c -o cellDynamics.so
mpicc $gcc_flags grid.c -o grid.so
mpicc $gcc_flags wind.c -o wind.so
mpicc $gcc_flags spatialDynamics.c -o spatialDynamics.so
mpicc $gcc_flags populationsAndParameters.c -o populationsAndParameters.so
mpicc $gcc_flags spatialDependence.c -o spatialDependence.so
