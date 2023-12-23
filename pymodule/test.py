import sys
import os
# import json
import numpy as np
import time
# from scipy.fft import fft
curdir = os.path.dirname(__file__)
cdir = os.path.abspath(os.path.join(curdir,'../mono_lags'))
sys.path.append(cdir)
import um21_lagmodel as um
# from matplotlib import pyplot
import example


class Time:
    def __init__(self):
        self._start = time.time()

    def timer(self):
        return time.time() - self._start

    def reset(self):
        self._start = time.time()


# class File:
#     def __init__(self, filename):
#         self._filename = filename

#         with open(os.path.join(os.getcwd(), filename), "r") as file:
#             self._data = json.load(file)

#             for function in self._data.values():
#                 for test in function.values():
#                     for put in test.values():
#                         for i, parameter in enumerate(put):
#                             put[i] = self.check(parameter)

#     def check(self, l):
#         if type(l) == list:
#             for i, el in enumerate(l):
#                 l[i] = self.check(el)
#         elif type(l) == dict and "real" in l and "imag" in l:
#             return complex(l["real"], l["imag"])
#         return l

#     @property
#     def data(self):
#         return self._data


def main(argv):
    # if len(argv) != 2:
    #     exit(f"Usage: {sys.argv[0]} <file> <threshold>")

    print(example.myadd(2, 3))
    print(example.myadd(3.3, 8.5))
    print(example.myadd("Hello ", "World"))

    for x in range(2, 6):
        n = pow(10, x)
        a = np.array(list(range(int(n))))
        a = np.divide(a, n/10)

        timer = Time()
        print(example.find_nearest(a, np.pi))
        print(timer.timer(), x, "C++")
        timer.reset()
        print(um.find_nearest(a, np.pi))
        print(timer.timer(), x, "python")

    # tresh = float(argv[1])

    # dbg = 0

    # CPP = File(argv[0])

    # functions = CPP.data.keys()

    # for function in functions:
    #     dbg += test_f(dbg, CPP, function, tresh)

    # # if len(um.test) != len(um.tset):
    # #     print(", ".join(i for i in um.tset if i not in um.test))
    # # else:
    # #     print("FINISHED!")

    # # print(dbg)
    # return dbg


# def test_f(dbg, CPP, function, thresh=0):
#     """Test the find_nearest function"""

#     test1_CPP = CPP.data[function]

#     assert type(function) == str

#     f = getattr(um, function)

#     for tset, test in test1_CPP.items():
#         for i in range(len(test["input"])):
#             if type(test["input"][i]) == list:
#                 test["input"][i] = np.array(test["input"][i])
#             if function == "calc_illumination_fracs" and \
#                     type(test["input"][i]) == str:
#                 test["input"][i] = getattr(um, test["input"][i])
            

#         test1 = f(*test["input"])

#         try:
#             assert len(test1) == len(test["output"])
#         except AssertionError:
#             assert len(test1) == len(test["output"][0])

#         if checker(test1, test["output"], thresh):
#             dbg += 1
#             print(f"ERROR in {tset} of {function}")
            
#             print("ERROR: ", end="")
#             try:
#                 print()
#                 test1[0][0]
#                 for i in range(len(test1)):
#                     print(f"python: {np.array(test1[i])}")
#                     print(f"C++:    {np.array(test['output'][i])}")
#             except IndexError:
#                 print(f"{test1} != {test['output']}")
#             exit()

#     return dbg


# def checker(data1, data2, thresh = 0):
#     try:
#         assert len(data1) == len(data2)
#     except AssertionError:
#         try:
#             assert len(data1) == len(data2[0])
#             data2 = data2[0]
#         except TypeError:
#             assert len(data1[0]) == len(data2)
#             data1 = data1[0]

#     for i in range(len(data1)):
#         if not (np.isscalar(data1[i]) and np.isscalar(data2[i])):
#             if checker(data1[i], data2[i], thresh):
#                 print("cheker fault")
#                 return True
#         elif abs(data1[i] - data2[i]) > thresh:
#             print(data1[i], data2[i], data1[i] - data2[i], ">", thresh)
#             print("not the same!")
#             return True

#     return False


# def geometry(r, parms):
#     return 0, 0*r, 0*r


if __name__ == "__main__":
    # setattr(um, "geometry", geometry)
    main(sys.argv[1:])