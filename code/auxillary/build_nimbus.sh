#!/bin/bash

cd auxillary; ./build_nimbus.sh ; cd ..

# get the python directory (eg /usr/include/python3.8/) by searching for Python.h on your system
# get the numpy include directory (eg /usr/lib/python3/dist-packages/numpy/core/include) using import numpy then numpy.get_include()
gcc_flags="-shared -fno-strict-aliasing -Wsign-compare -Wunreachable-code -DNDEBUG -g -O3 -flto -Wl,--rpath -Wl,${LD_RUN_PATH} -Wall -I/usr/include/python3.8/ -I/usr/lib/python3/dist-packages/numpy/core/include -fPIC -march=native"

gcc $gcc_flags csvparser.c -o csvparser.so

cython -3 -a *.pyx

gcc $gcc_flags cellDynamics.c -o cellDynamics.so
gcc $gcc_flags grid.c -o grid.so
gcc $gcc_flags wind.c -o wind.so
gcc $gcc_flags spatialDynamics.c -o spatialDynamics.so
gcc $gcc_flags populationsAndParameters.c -o populationsAndParameters.so
gcc $gcc_flags spatialDependence.c -o spatialDependence.so

