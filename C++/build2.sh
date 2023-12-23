#!/bin/bash

c++ -O3 -Wall -fno-inline -shared -std=c++17 -fPIC $(python -m pybind11 --includes) pybind*.cpp -o lagmodel$(python3-config --extension-suffix) -Wno-sign-compare
c++ -O3 -Wall -fno-inline -std=c++17 -g copy_check.cpp -o copy_check $(python -m pybind11 --includes) -Wno-sign-compare  -lpython3.8
# valgrind  --vgdb=yes --vgdb-error=0 ./copy_check
valgrind ./copy_check --leak-check=full