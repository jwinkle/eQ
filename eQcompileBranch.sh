#!/bin/bash

mkdir -p build
rm -rf ./build/*

cd build
cmake ..
make

# ./eQ test
/usr/bin/mpirun -n 2 -display-allocation  ./eQ test
