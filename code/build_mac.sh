#!/bin/bash

cd auxillary; ./build_mac.sh ; cd ..

cythonize -3 -f -i -a *.pyx
