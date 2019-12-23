#!/bin/bash

gcc_flags="-shared -fno-strict-aliasing -Wsign-compare -Wunreachable-code -DNDEBUG -g -O3 -Wall -fPIC"

mpicc -g -O3 -fno-strict-aliasing -Wsign-compare -Wunreachable-code -DNDEBUG -Wall  -c ../csvparser.c
mpicc -dynamiclib csvparser.o -install_name @rpath/t/libcsvparser.dylib -o libcsvparser.dylib
mkdir -p t ; mv libcsvparser.dylib t/
mpicc -Wall -lcsvparser -L`pwd`/t -Xlinker -rpath -Xlinker `pwd` ab_convert.c -o ab_convert


