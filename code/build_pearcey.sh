#!/bin/bash

cythonize -3 -a *.pyx

gcc_flags="-shared -fno-strict-aliasing -Wsign-compare -Wunreachable-code -DNDEBUG -g -O3 -Wall -I/apps/python/3.7.2/include/python3.7m -fPIC"
mpicc $gcc_flags grid.c -o grid.so
mpicc $gcc_flags wind.c -o wind.so
mpicc $gcc_flags cell.c -o cell.so
mpicc $gcc_flags spatial.c -o spatial.so
mpicc $gcc_flags populations.c -o populations.so
mpicc $gcc_flags spatialDependence.c -o spatialDependence.so
