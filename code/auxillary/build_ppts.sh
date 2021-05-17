#!/bin/bash

module load openmpi/4.1.1-ofed51
module load git
module load python/3.9.4

gcc_flags="-shared -fno-strict-aliasing -Wsign-compare -Wunreachable-code -DNDEBUG -g -O3 -flto -Wl,--rpath -Wl,${LD_RUN_PATH} -Wall -fPIC"

mpicc $gcc_flags ../csvparser.c -o libcsvparser.so
mpicc -Wall -L. -lcsvparser -Wl,-rpath=${PWD} ab_convert.c -o ab_convert


