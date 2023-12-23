#!/bin/bash

c++ -O3 -Wall -shared -std=c++17 -fPIC $(python -m pybind11 --includes) *.cpp -o example$(python3-config --extension-suffix)