import sys
import os
import json
import numpy as np
from scipy.fft import fft
sys.path.append(os.path.abspath("../"))
from mono_lags import um21_lagmodel as um
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

def main(argv):
    if len(argv) != 0:
        exit(f"Usage: {sys.argv[0]}")

    N = int(1E3)
    N = 400

    print("Difference: Postive = python slower, Negative = C++ slower")
    print("Ratio: bigger than 1 = python slower, smaller than 1 = C++ slower")

    t = Time()
    a = np.linspace(0, 1, N)
    b = a * 10
    create = t.timer(True)

    print("geometry")
    t.timer(True)

    # um.calc_illumination_fracs(a, b, getattr(um, "geometry"), [0.37])
    # py = t.timer(True)
    # print(f"python: {py:.2f}")
    # lm.calc_illumination_fracs(geo.Geometry(a, b, .37))
    # cpp = t.timer(True)
    # print(f"C++: {cpp:.2f}")

    # print(f"Difference: {py - cpp + create:.2f}; Ratio: {(py+create)/cpp:.2f}")
        
    # print("an_shere")
    # t.timer(True)

    # um.calc_illumination_fracs(a, b, getattr(um, "an_sphere"), [0.37])
    # py = t.timer(True)
    # print(f"python: {py:.2f}")
    # lm.calc_illumination_fracs(geo.AN_Sphere(a, b, .37))
    # cpp = t.timer(True)
    # print(f"C++: {cpp:.2f}")

    # print(f"Difference: {py - cpp + create:.2f}; Ratio: {(py+create)/cpp:.2f}")
    
    # parms = [.13, .37, .7, 6.7/2, 6.7/7, 6.7/11, 6.7]

    # print("bknpow_emiss")
    # t.timer(True)

    # um.calc_illumination_fracs(a, b, getattr(um, "bknpow_emiss"), parms)
    # py = t.timer(True)
    # print(f"python: {py:.2f}")
    # lm.calc_illumination_fracs(geo.BKNpow_Emiss(a, b, *parms))
    # cpp = t.timer(True)
    # print(f"C++: {cpp:.2f}")

    # print(f"Difference: {py - cpp + create:.2f}; Ratio: {(py+create)/cpp:.2f}")
    
    # parms = [.13, .7, .37, 6.7/2, 6.7/7, 6.7/11, 6.7]

    # print("bknpow_emiss")
    # t.timer(True)

    # um.calc_illumination_fracs(a, b, getattr(um, "bknpow_emiss"), parms)
    # py = t.timer(True)
    # print(f"python: {py:.2f}")
    # lm.calc_illumination_fracs(geo.BKNpow_Emiss(a, b, *parms))
    # cpp = t.timer(True)
    # print(f"C++: {cpp:.2f}")

    # print(f"Difference: {py - cpp + create:.2f}; Ratio: {(py+create)/cpp:.2f}")
    
    parms = [.37, .37*3, 400, 100]

    print(f"cylinder ({'x'.join(map(str, parms[-2:]))})")
    t.timer(True)

    um.calc_illumination_fracs(a, b, getattr(um, "cylinder"), parms)
    py = t.timer(True)
    print(f"python: {py:.2f}")
    lm.calc_illumination_fracs(geo.Cylinder(a, b, *parms))
    cpp = t.timer(True)
    print(f"C++: {cpp:.2f}")

    print(f"Difference: {py - cpp + create:.2f}; Ratio: {(py+create)/cpp:.2f}")

    # parms = [.37, .37*3, 21, 21]

    # print(f"cylinder ({'x'.join(map(str, parms[-2:]))})")
    # t.timer(True)

    # um.calc_illumination_fracs(a, b, getattr(um, "cylinder"), parms)
    # py = t.timer(True)
    # print(f"python: {py:.2f}")
    # lm.calc_illumination_fracs(geo.Cylinder(a, b, *parms))
    # cpp = t.timer(True)
    # print(f"C++: {cpp:.2f}")

    # print(f"Difference: {py - cpp + create:.2f}; Ratio: {(py+create)/cpp:.2f}")

    # parms = [.37, .37*3, 47, 47]

    # print(f"cylinder ({'x'.join(map(str, parms[-2:]))})")
    # t.timer(True)

    # um.calc_illumination_fracs(a, b, getattr(um, "cylinder"), parms)
    # py = t.timer(True)
    # print(f"python: {py:.2f}")
    # lm.calc_illumination_fracs(geo.Cylinder(a, b, *parms))
    # cpp = t.timer(True)
    # print(f"C++: {cpp:.2f}")

    # print(f"Difference: {py - cpp + create:.2f}; Ratio: {(py+create)/cpp:.2f}")
    
    # parms = [.13, .37, .7, 6, 6]

    # print(f"inv_cone ({'x'.join(map(str, parms[-2:]))})")
    # t.timer(True)

    # um.calc_illumination_fracs(a, b, getattr(um, "inv_cone"), parms)
    # py = t.timer(True)
    # print(f"python: {py:.2f}")
    # lm.calc_illumination_fracs(geo.Inv_Cone(a, b, *parms))
    # cpp = t.timer(True)
    # print(f"C++: {cpp:.2f}")

    # print(f"Difference: {py - cpp + create:.2f}; Ratio: {(py+create)/cpp:.2f}")

    # parms = [.13, .37, .7, 21, 21]

    # print(f"inv_cone ({'x'.join(map(str, parms[-2:]))})")
    # t.timer(True)

    # um.calc_illumination_fracs(a, b, getattr(um, "inv_cone"), parms)
    # py = t.timer(True)
    # print(f"python: {py:.2f}")
    # lm.calc_illumination_fracs(geo.Inv_Cone(a, b, *parms))
    # cpp = t.timer(True)
    # print(f"C++: {cpp:.2f}")

    # print(f"Difference: {py - cpp + create:.2f}; Ratio: {(py+create)/cpp:.2f}")

    # parms = [.13, .37, .7, 47, 47]

    # print(f"inv_cone ({'x'.join(map(str, parms[-2:]))})")
    # t.timer(True)

    # um.calc_illumination_fracs(a, b, getattr(um, "inv_cone"), parms)
    # py = t.timer(True)
    # print(f"python: {py:.2f}")
    # lm.calc_illumination_fracs(geo.Inv_Cone(a, b, *parms))
    # cpp = t.timer(True)
    # print(f"C++: {cpp:.2f}")

    # print(f"Difference: {py - cpp + create:.2f}; Ratio: {(py+create)/cpp:.2f}")
    
    # parms = [.37, 6, 6]

    # print(f"sphere ({'x'.join(map(str, parms[-2:]))})")
    # t.timer(True)

    # um.calc_illumination_fracs(a, b, getattr(um, "sphere"), parms)
    # py = t.timer(True)
    # print(f"python: {py:.2f}")
    # lm.calc_illumination_fracs(geo.Sphere(a, b, *parms))
    # cpp = t.timer(True)
    # print(f"C++: {cpp:.2f}")

    # print(f"Difference: {py - cpp + create:.2f}; Ratio: {(py+create)/cpp:.2f}")

    # parms = [.37, 21, 21]

    # print(f"sphere ({'x'.join(map(str, parms[-2:]))})")
    # t.timer(True)

    # um.calc_illumination_fracs(a, b, getattr(um, "sphere"), parms)
    # py = t.timer(True)
    # print(f"python: {py:.2f}")
    # lm.calc_illumination_fracs(geo.Sphere(a, b, *parms))
    # cpp = t.timer(True)
    # print(f"C++: {cpp:.2f}")

    # print(f"Difference: {py - cpp + create:.2f}; Ratio: {(py+create)/cpp:.2f}")

    # parms = [.37, 47, 47]

    # print(f"sphere ({'x'.join(map(str, parms[-2:]))})")
    # t.timer(True)

    # um.calc_illumination_fracs(a, b, getattr(um, "sphere"), parms)
    # py = t.timer(True)
    # print(f"python: {py:.2f}")
    # lm.calc_illumination_fracs(geo.Sphere(a, b, *parms))
    # cpp = t.timer(True)
    # print(f"C++: {cpp:.2f}")

    # print(f"Difference: {py - cpp + create:.2f}; Ratio: {(py+create)/cpp:.2f}")


def geometry(r, parms):
    return 0, 0*r, 0*r


if __name__ == "__main__":
    setattr(um, "geometry", geometry)
    main(sys.argv[1:])


# nr = 100
# nr_area = 100
# meshgrid = 1000 x 1000

# nr = 400
# nphi = 100
# nz = 400