#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include "functions.h"

namespace py = pybind11;

PYBIND11_MODULE(lagmodel, m) {
    m.doc() = "Module with functions for the Lag model";
    m.def("find_nearest", &find_nearest<double>,
        "find nearest value in array and report it's index;\n\n"
        "find_nearest(array, value) -> (value, index)");
    m.def("calc_dispfrac", &calc_dispfrac<double>,
        "calculate the disspation fraction;\n\n"
        "calc_dispfrac(rad, rad_area, rin, rcor, seedff_norm, seedff_ind, "
        "heatff_norm, heatff_ind) -> (dispfrac, seef_frac_flow, "
        "heat_frac_flow)");
    m.def("calc_illumination_fracs", &calc_illumination_fracs<double>,
        "calculate the illumination fractions;\n\n"
        "calc_illumination_fracs(geomod) -> (omega_cor, frad_disktocor, "
        "frad_cortodisk)");
    m.def("calc_timing_params", &calc_timing_params<double>,
        "calculate the timing parameters;\n\n"
        "calc_timing_params(rad, i_rsigmax, rcor, i_rcor, t_scale, "
        "disk_tau_par, cor_tau_par, lor_model, lor_par) -> (tau, lfreq, "
        "q_list, rms)");
    m.def("calc_propagation_params", &calc_propagation_params<double>,
        "calculate the propagation parameters;\n\n"
        "calc_propagation_params(rad, rad_edge, rcor, disk_prop_par, "
        "cor_prop_par) -> (deltau)");

    auto internal = m.def_submodule("internal_class");
    // internal.def()

    auto geo = m.def_submodule("geometries", "Geometry classes for the lagmodel"
        " module");

    typedef Geometry<list, double> geo_base;

    py::class_<geo_base>(geo, "Geometry")
        .def(py::init<list<double>, list<double>, double>())
        .def("get_rcor", &geo_base::get_rcor)
        .def("get_omega_cor", &geo_base::get_omega_cor)
        .def("get_frad_disktocor", &geo_base::get_frad_disktocor)
        .def("get_frad_cortodisk", &geo_base::get_frad_cortodisk)
        .def("get_omega_cor_area", &geo_base::get_omega_cor_area)
        .def("get_frad_disktocor_area", &geo_base::get_frad_disktocor_area)
        .def("get_frad_cortodisk_area", &geo_base::get_frad_cortodisk_area);
}