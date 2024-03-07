#!/bin/bash

cd auxillary; ./build_mac.sh ; cd ..

cythonize -3 -f -i -a --directive=profile=True --directive=linetrace=True --option=CYTHON_TRACE=1 *.pyx
