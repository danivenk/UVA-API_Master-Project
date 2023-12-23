import sys
import os
import json
import numpy as np
from scipy.fft import fft
# from mono_lags import um_lagmodel as um
from matplotlib import pyplot
import time

import lagmodel as lm
from lagmodel import geometries as geo, internal_class as ic

class Time:
    def __init__(self):
        self._start = time.time()

    def timer(self, reset=False):
        finish = time.time() - self._start

        if reset:
            self._start = time.time()

        return finish

class File:
    def __init__(self, filename):
        self._filename = filename

        with open(os.path.join(os.getcwd(), filename), "r") as file:
            self._data = json.load(file)

            for function in self._data.values():
                for test in function.values():
                    for put in test.values():
                        for i, parameter in enumerate(put):
                            put[i] = self.check(parameter)

    def check(self, l):
        if type(l) == list:
            l = np.array(l)
            for i, el in enumerate(l):
                l[i] = self.check(el)
        elif type(l) == dict and "real" in l and "imag" in l:
            return complex(l["real"], l["imag"])
        return l

    @property
    def data(self):
        return self._data

def main(argv):
    if len(argv) != 1:
        exit(f"Usage: {sys.argv[0]} <filename>")

    # os.system("./build.sh")

    CPP = File(argv[0]).data
    func = list(CPP.keys())

    cname = {"geometry": "Geometry", "an_sphere": "AN_Sphere",
        "bknpow_emiss": "BKNpow_Emiss", "cylinder": "Cylinder",
        "inv_cone": "Inv_Cone", "sphere": "Sphere",
        "piecewise_emiss": "Piecewise_Emiss"}

    print("Difference: Postive = python slower, Negative = C++ slower")
    print("Ratio: bigger than 1 = python slower, smaller than 1 = C++ slower")
    print("; ".join(f"{i}: {item}" for i, item in enumerate(func)))

    t = Time()

    cf = func[13]

    test = CPP[cf]["test_1"]

    print(5*"-", cf)

    i = 0

    # print(test["input"])
    print(test["output"][0])
    # a = getattr(lm, cf)(*test["input"], True)
    # a = getattr(lm, cf)(getattr(geo, cname[test["input"][2]])(*test["input"][0:2], *tuple(*test["input"][3:])))
    a = getattr(lm, "lorentz_q")(*test["input"])
    print(a)


def geometry(r, parms):
    return 0, 0*r, 0*r


if __name__ == "__main__":
    # setattr(um, "geometry", geometry)
    main(sys.argv[1:])