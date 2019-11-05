#!/bin/bash

cythonize -3 -a *.pyx

gcc_flags="-shared -fno-strict-aliasing -Wsign-compare -Wunreachable-code -DNDEBUG -g -O3 -Wall -I/apps/python/3.7.2/include/python3.7m -fPIC"
mpicc $gcc_flags grid.c -o grid.so

