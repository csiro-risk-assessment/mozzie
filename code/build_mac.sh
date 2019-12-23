#!/bin/bash

cd auxillary; ./build_mac.sh ; cd ..

cythonize -3 -i -a *.pyx
