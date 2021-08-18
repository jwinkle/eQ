#!/bin/bash

# COMPILE FENICS .UFL FILES INTO .H USING FFC
cd fenics
./ufl_ffc.sh

cd ..

#MAKE (CLEAN) BUILD DIRECTORY AND CD INTO; RUN CMAKE (CMakeLists.txt file is in root dir); build project
mkdir -p build
cd build
rm -rf [Cc][Mm]ake*

cmake ..
make

# ./eQ test using 3 MPI nodes:
/usr/bin/mpirun -n 3  ./eQ test
