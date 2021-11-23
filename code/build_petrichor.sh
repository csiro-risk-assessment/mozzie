#!/bin/bash

module load openmpi/4.1.1-ofed51
module load gcc
module load git
module load python/3.9.4


cd auxillary; ./build_ppts.sh ; cd ..

# get the numpy include directory using import numpy then numpy.get_include()
gcc_flags="-shared -fno-strict-aliasing -Wsign-compare -Wunreachable-code -DNDEBUG -g -O3 -flto -Wl,--rpath -Wl,${LD_RUN_PATH} -Wall -I/apps/python/3.9.4/include/python3.9 -I/apps/python/3.9.4/lib/python3.9/site-packages/numpy-1.20.2-py3.9-linux-x86_64.egg/numpy/core/include -fPIC -march=native"

mpicc $gcc_flags csvparser.c -o csvparser.so

cythonize -3 -a *.pyx

mpicc $gcc_flags cellDynamics.c -o cellDynamics.so
mpicc $gcc_flags grid.c -o grid.so
mpicc $gcc_flags wind.c -o wind.so
mpicc $gcc_flags spatialDynamics.c -o spatialDynamics.so
mpicc $gcc_flags populationsAndParameters.c -o populationsAndParameters.so
mpicc $gcc_flags spatialDependence.c -o spatialDependence.so
