#!/bin/bash

cd davidSquarer && ./machconfig && make && cd ..
cd ilkkaSquarer && make && cd ..
cd exact && f2py -c -m coulcc coulcc36-f90.f90 && cd ..
cd ewald && [[ -d builds ]] || mkdir builds && cd build && cmake .. && make && cd .. && cd ..
