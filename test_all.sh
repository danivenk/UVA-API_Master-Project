#!/bin/bash
cd C++; make clean; make; valgrind ./main
cd ../; python test.py C++/cpp_out.txt 7.5E-15