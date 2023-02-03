#!/bin/bash

c++ -O3 -Wall -shared -std=c++17 -fPIC $(python -m pybind11 --includes) pybind*.cpp -o lagmodel$(python3-config --extension-suffix) -Wno-sign-compare
# c++ -O3 -Wall -std=c++17 test.cpp -o test $(python -m pybind11 --includes) -Wno-sign-compare
./test