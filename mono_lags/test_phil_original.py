import os
import sys
import numpy as np
sys.path.append(os.path.abspath("../C++/"))
import lagmodel as lm
from lagmodel import geometries as geo
import um21_lagmodel_new as uml

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

    string = argv

    params = {}
    parms = string.split("), (")
    parms[0] = parms[0][13:]
    parms[-1] = parms[-1][:-4]

    for i, parm in enumerate(parms):
        parm = Parameter(parm)
        params[parm.name] = parm

    rin = params['rin'].value
    rcor = params['rcor'].value
    hcor = params['hcor'].value
    cone_angle = params['cone_angle'].value
    disk_tau_norm = params['disk_tau_norm'].value   
    disk_tau_ind = params['disk_tau_ind'].value

    freq_rsigmax=0.01
    t_scale=5e-5
    rout = 1000
    nrad = 100
    nphi = 1000
    nz = 1000

    rtop = rcor + hcor*np.tan(np.pi*(np.array(cone_angle)/180))
    geopar = [rcor,hcor,rtop,nphi,nz]

    rsigmax = ((1/(freq_rsigmax*t_scale))/(disk_tau_norm*rcor**(-1.*disk_tau_ind)))**(1/(1.5+disk_tau_ind))
    
    rout = max([2*rsigmax,rout])
        
    #First set up radial grid
    rad_edge1 = np.logspace(np.log10(rin), np.log10(rsigmax), nrad+1)  # set up radial bin edges
    rad_edge2 = np.logspace(np.log10(rsigmax), np.log10(rout), nrad+1)
    rad_edge = np.append(rad_edge1,rad_edge2[1:])
    nrad = nrad*2
    
    rad = np.sqrt(rad_edge[1:]*rad_edge[:-1])  # Define radial bin centres as geometric mean of bin edges
    rad_area = np.pi*(np.square(rad_edge[1:])-np.square(rad_edge[:-1]))
    # Now reset coronal outer radius (i.e. disk inner radius) to be at the nearest radial bin edges
    rcor, i_rcor = uml.find_nearest(rad_edge,rcor)

    geometry = geo.Inv_Cone(rad, rad_area, *geopar)
    omega_cor, disktocor_frac, cortodisk_frac = lm.calc_illumination_fracs(geometry)
    omega_cor1, disktocor_frac1, cortodisk_frac1 = uml.calc_illumination_fracs(rad,rad_area,uml.inv_cone,geopar)

    return

if __name__ == "__main__":
    if len(sys.argv) < 1:
        print(f"Usage: {sys.argv[0]} <filename>")
        exit(sys.argv)
    for argv in sys.argv[1:]:
        main(argv)