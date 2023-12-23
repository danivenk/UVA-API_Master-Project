import os
import sys
import numpy as np

class Parameter:
    def __init__ (self, string):
        vals = string.split(", ")
        self._name = vals[0][1:-1]
        if "fixed" in vals[2]:
            self._value = float(vals[2].split("=")[1].split(" (")[0])
        else:
            self._value = float(vals[2].split("=")[1])
        self._bounds = vals[3][:-1].split("=")[1].split(":")
        self._bounds[0] = self._bounds[0][1:]
        self._bounds[1] = self._bounds[1][:-1]
        try:
            assert len(self._bounds) == 2
        except:
            print(self._bounds)
            raise AssertionError("length of bounds was not 2")
        for i, bound in enumerate(self._bounds):
            self._bounds[i] = float(bound)

    @property
    def name(self):
        return self._name
    
    @property
    def value(self):
        return self._value
    
    @property
    def bounds(self):
        return self._bounds

    @property
    def bounds_size(self):
        return np.diff(self._bounds)[0]

    def __sub__(self, other):
        if self == other:
            return None
        if self._name != self._name:
            return AssertionError("not the same type of parameter")
        return self._value - other._value
    
    def __truediv__(self, other):
        if self == other:
            return None
        if self._name != self._name:
            return AssertionError("not the same type of parameter")
        return self._value/other._value

def main(argv):
    print("Current file", argv)

    # phil >0; dani <0
    chi = []
    params = []

    with open(argv) as file:
        for line in file:
            if "Parameters" in line:
                parms = line.split("), (")
                parms[0] = parms[0][13:]
                parms[-1] = parms[-1][:-4]

                for i, parm in enumerate(parms):
                    parm = Parameter(parm)
                    try:
                        params[i].append(parm)
                    except IndexError:
                        params.append([parm])
            else:
                try:
                    chi.append(float(line))
                except ValueError:
                    pass

    chid = np.diff(chi)

    for parm in params:
        name = parm[0]._name
        parmd = np.diff(parm)
        if np.sum(parmd) == 0:
            continue
        print("\t" + 10*"=", name, 10*"=")
        for i in range(len(chid)):
            p = (1 - parm[i+1]/parm[i])*100
            pnew = (parmd[i]/parm[0].bounds_size)*100
            try:
                c = (1 - chi[i+1]/chi[i])*100
            except ZeroDivisionError:
                c = (1 - chi[i+1]/(chi[i]+1E-100))*100
            if abs(c) > 50:
                print(f"\t\t{i:<4}: {pnew:.1f} # {chid[i]:.1f} ({c:.1f}%)")
                # print(f"\t\t{i:<4}: {parmd[i]:.1f} ({p:.1f}%) # "
                #       f"{chid[i]:.1f} ({c:.1f}%)")
                print(f"\t\t{parm[i-1]._value} - {chi[i-1]}")
                print(f"\t\t{parm[i]._value} - {chi[i]}")
                print(f"\t\t{parm[i+1]._value} - {chi[i+1]}")

    return

if __name__ == "__main__":
    
    if len(sys.argv) <= 1:
        print(f"Usage: {sys.argv[0]} <filename>")
        exit(sys.argv)
    for arg in sys.argv[1:]:
        main(arg)