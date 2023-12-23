import sys
import os
import json
import numpy as np
from scipy.fft import fft
from mono_lags import um21_lagmodel as um
from matplotlib import pyplot
import time

sys.path.append(os.path.abspath("C++/"))
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

    CPP = File(argv[0])

    cname = {"geometry": "Geometry", "an_sphere": "AN_Sphere",
        "bknpow_emiss": "BKNpow_Emiss", "cylinder": "Cylinder",
        "inv_cone": "Inv_Cone", "sphere": "Sphere"}

    print("Difference: Postive = python slower, Negative = C++ slower")
    print("Ratio: bigger than 1 = python slower, smaller than 1 = C++ slower")
    print(CPP.data.keys())

    t = Time()

    for f in CPP.data.keys():
        if f == "calc_illumination_fracs":
            # print(CPP.data[f])
            for test, values in CPP.data[f].items():
                values["inputc"] = values["input"]
                parms = values["input"][3]
                geom = getattr(geo, cname[values["input"][2]])
                values["inputc"] = geom(*values["input"][0:2], *parms)
            # print(CPP.data[f])
        print(f)
        for test, values in CPP.data[f].items():
            parms = values["input"]
            if f == "calc_illumination_fracs":
                # continue
                print(test, parms[2])
                parmsc = values["inputc"]
                parms[2] = getattr(um, parms[2])
            else:
                print(test)
            if f == "lorentz_qold":
                func_cpp = getattr(lm, "lorentz_q")
            else:
                func_cpp = getattr(lm, f)
            func_py = getattr(um, f)
            t.timer(True)
            func_py(*parms)
            py = t.timer(True)
            if f == "calc_illumination_fracs":
                func_cpp(parmsc)
            else:
                func_cpp(*parms)
            cpp = t.timer(True)
            print(f"Difference: {py - cpp:.2E}; Ratio: "
                  f"{py/cpp:.2f}")


def geometry(r, parms):
    return 0, 0*r, 0*r


if __name__ == "__main__":
    setattr(um, "geometry", geometry)
    main(sys.argv[1:])