#!/bin/bash

cd davidSquarer && ./machconfig && make && cd -
cd ilkkaSquarer && make && cd -
cd exact && f2py -c -m coulcc coulcc36-f90.f90 && cd -
