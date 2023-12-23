import sys
import os
import json
import numpy as np
from scipy.fft import fft
# import mono_lags.um21_lagmodel as um
from matplotlib import pyplot

from lagmodel.internal_class import Array, Nested_Array
import time


def main(argv):
    N = int(argv[0])
    print(f"No of items: {N*N}")
    start = time.time()
    Nested_Array(N,N).get_element(34, 87)
    end = time.time()
    f1 = end-start
    print(f"Nested_Array: {end-start:.2f}s")
    start = time.time()
    Array(N,N).get_element(35, 87)
    end = time.time()
    f2 = end-start
    print(f"Array: {end-start:.2f}s")

    print(f"Difference: {abs(f1-f2):.2f}; Ratio: {f1/f2:.2f}")


if __name__ == "__main__":
    main(sys.argv[1:])

# try geometry with 1M compare to python