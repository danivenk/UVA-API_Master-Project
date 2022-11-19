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
    m.def("calc_radial_time_response", &calc_radial_time_response<double>,
        "calculate the radial time response;\n\n"
        "calc_radial_time_repsonse(rad, rcor, disp_frac, disktocor_frac, "
        "cortodisk_frac, seed_frac_flow, heat_frac_flow) -> (ldisk_disp, "
        "lseed_disp, lheat, ldisk_rev, lseed_rev)");
    m.def("calc_irfs_mono", &calc_irfs_mono<double, Array<list, double>>,
        "caclulate the irfs mono;\n\n"
        "calc_irfs_mono(gama_par, e_seed, energy, ldisk_disp, lseed_disp, "
        "lheat, ldisk_rev, lseed_rev) -> (gamma_mean, gamma_irf, flux_irf, "
        "disk_irf, seed_irf)");
    // m.def("calc_irfs_mono", &calc_irfs_mono<double, Nested_Array<list, double>>);
    m.def("linear_rebin_irf", &linear_rebin_irf<double>,
        "rebin to linear scale;\n\n"
        "linear_rebin_irf(dt, i_rsigmax, irf_nbins, irf_binedgefrac, "
        "input_irf, deltau_scale, nirf) -> (rebinned_irf, flux_outer)");
    m.def("calc_cross_psd", &calc_cross_psd<double>,
        "calculate the cross power spectra\n\n"
        "calc_cross_psd(freq, dt, ci_irf, ref_irf, irf_nbis, deltau_scale, "
        "i_rsigmax, f_pk, q, rms) -> (wt_cross_spec, wt_pow_spec_ci, "
        "wt_pow_spec_ref, mod_sig_psd");
    m.def("lorentz_q", &lorentz_q<double, Array<list, double>>,
        "calculates the lorentzian function;\n\n"
        "lorentz_q(f, f_pk, q, rms) -> (lorentz)");
    // m.def("lorentz_q", &lorentz_q<double, Nested_Array<list, double>>);
    m.def("lorentz_q", &lorentz_q<double>);
    m.def("calculate_stprod_mono", &calculate_stprod_mono<double,
        Array<list, double>, Array<list, int>>,
        "calculates mono-energetic spectral-timing products;\n\n"
        "calculate_stprod_mono(nirf_mult, energy, encomb, flux_irf, disk_irf, "
        "gamma_irf, deltau, min_deltau_frac, i_rsigmax, lfreq, q, rms, "
        "t_scale) -> (freq, phlag, tlag, psd_ci, psd_ref, "
        "mod_sig_psd,irf_nbins, irf_binedgefrac, deltau_scale, dt, nirf, "
        "ci_irf, ci_mean, ci_outer)");

    auto internal = m.def_submodule("internal_class");
    
    typedef Matrix<list, double> matrix_0;

    py::class_<matrix_0>(internal, "Matrix")
        .def(py::init<>())
        .def("get_size", &matrix_0::get_size);

    typedef Nested_Array<list, double> matrix_1;

    py::class_<matrix_1, matrix_0>(internal, "Nested_Array")
        .def(py::init<list<double>, list<double>>())
        .def(py::init<matrix_1>())
        .def(py::init<int, int>())
        .def("get_element", py::overload_cast<int, int>(&matrix_1::get_element))
        .def("get_element", py::overload_cast<int, int>(&matrix_1::get_element,
            py::const_))
        .def("get_matrix", &matrix_1::get_matrix);

    typedef Array<list, double> matrix_2;

    py::class_<matrix_2, matrix_0>(internal, "Array")
        .def(py::init<list<double>, list<double>>())
        .def(py::init<matrix_2>())
        .def(py::init<int, int>())
        .def("get_element", py::overload_cast<int, int>(&matrix_2::get_element))
        .def("get_element", py::overload_cast<int, int>(&matrix_2::get_element,
            py::const_))
        .def("get_matrix", &matrix_2::get_matrix);

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