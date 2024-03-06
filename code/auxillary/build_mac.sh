#!/bin/sh

cc -g -O3 -fno-strict-aliasing -Wsign-compare -Wunreachable-code -DNDEBUG -Wall  -c ../csvparser.c
cc -dynamiclib csvparser.o -install_name @rpath/t/libcsvparser.dylib -o libcsvparser.dylib
mkdir -p t ; mv libcsvparser.dylib t/
cc -Wall -lcsvparser -L`pwd`/t -Xlinker -rpath -Xlinker `pwd` ab_convert.c -o ab_convert


