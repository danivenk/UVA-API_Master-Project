#!/bin/bash

c++ -O3 -Wall -g -shared -std=c++17 -fPIC $(python -m pybind11 --includes) timing_pybind*.cpp -o timing_lagmodel$(python3-config --extension-suffix) -Wno-sign-compare