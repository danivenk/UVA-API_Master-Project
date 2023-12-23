#!/bin/bash
g++ -O3 -Wall -Werror -shared -std=c++17 -fPIC *.cpp -o libcppmult.so
g++ -O3 -Wall -Werror -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` -I /usr/include/python3.7 -I . *.cpp -o library.so `python3.7-config --extension-suffix` "
    "-L. -lcppmult -Wl,-rpath,.".format(cpp_name, extension_name)