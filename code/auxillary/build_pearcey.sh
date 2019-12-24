#!/bin/bash

module load openmpi/3.1.4-ofed45-gcc
module load git
module load python/3.7.2

gcc_flags="-shared -fno-strict-aliasing -Wsign-compare -Wunreachable-code -DNDEBUG -g -O3 -Wall -fPIC"

mpicc $gcc_flags ../csvparser.c -o libcsvparser.so
mpicc -Wall -L. -lcsvparser -Wl,-rpath=${PWD} ab_convert.c -o ab_convert


