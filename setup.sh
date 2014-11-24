#!/bin/bash

git submodule init
git submodule update

cd davidSquarer && ./machconfig && make && cd ..
cd ilkkaSquarer && make && cd ..
cd exact && f2py -c -m coulcc coulcc36-f90.f90 && cd ..
cd ewald && [[ -d build ]] || mkdir build && cd build && cmake .. && make && cd .. && cd ..
