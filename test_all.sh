#!/bin/bash
cd C++; make clean; make; ./main
cd ../; python test.py C++/cpp_out.txt