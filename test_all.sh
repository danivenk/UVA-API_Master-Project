#!/bin/bash
cd C++; make clean; make; valgrind ./main
cd ../; python test.py C++/cpp_out.txt 1E-14