import sys
import os
import json
import numpy as np
from scipy.fft import fft
import mono_lags.um21_lagmodel as um
from matplotlib import pyplot

class File:
    def __init__(self, filename):
        self._filename = filename

        with open(os.path.join(os.getcwd(), filename), "r") as file:
            self._data = json.load(file)

    @property
    def data(self):
        return self._data


def main(argv):
    if len(argv) != 2:
        exit(f"Usage: {sys.argv[0]} <file> <threshold>")

    tresh = float(argv[1])

    dbg = 0

    CPP = File(argv[0])

    functions = CPP.data.keys()

    for function in functions:
        dbg += test_f(dbg, CPP, function, tresh)

    if len(um.test) != 0:
        print(", ".join(i for i in um.tset if i not in um.test))
    else:
        print("FINISHED!")

    # print(dbg)
    return dbg


def test_f(dbg, CPP, function, thresh=0):
    """Test the find_nearest function"""

    test1_CPP = CPP.data[function]

    assert type(function) == str

    f = getattr(um, function)

    for tset, test in test1_CPP.items():
        for i in range(len(test["input"])):
            if type(test["input"][i]) == list:
                test["input"][i] = np.array(test["input"][i])
            if function == "calc_illumination_fracs" and \
                    type(test["input"][i]) == str:
                test["input"][i] = getattr(um, test["input"][i])
            

        test1 = f(*test["input"])

        assert len(test1) == len(test["output"])

        if checker(test1, test["output"], thresh):
            dbg += 1
            print(f"ERROR in {tset}")
            
            print("ERROR: ")
            try:
                test1[0][0]
                for i in range(len(test1)):
                    print(f"python: {np.array(test1[i])}")
                    print(f"C++:    {np.array(test['output'][i])}")
            except IndexError:
                print(f"ERROR: {test1} != {test['output']}")

    return dbg


def checker(data1, data2, thresh = 0):
    assert len(data1) == len(data2)

    for i in range(len(data1)):
        if not (np.isscalar(data1[i]) and np.isscalar(data2[i])):
            if checker(data1[i], data2[i], thresh):
                print("cheker fault")
                return True
        elif abs(data1[i] - data2[i]) > thresh:
            print(data1[i], data2[i], data1[i] - data2[i], ">", thresh)
            print("not the same!")
            return True

    return False


def geometry(r, parms):
    return 0, 0*r, 0*r


if __name__ == "__main__":
    setattr(um, "geometry", geometry)
    main(sys.argv[1:])