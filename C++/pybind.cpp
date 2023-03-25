#include <variant>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include "functions.h"
// #include "pybind_cast.h"

namespace py = pybind11;
using namespace py::literals;

template<class T, class U>
py::array_t<T> as_array(U &m) {
    auto [row, col] = m.get_size();
    auto ptr_val = m.get_matrix().data();
    auto data_copy = std::make_unique<T[]>(row * col);
    std::memcpy(data_copy.get(), ptr_val, row * col * sizeof(T));
    return py::array_t<T>({row, col}, {sizeof(T) * col, sizeof(T)},
        data_copy.get());
}

template<class T, class U>
T as_matrix(py::array_t<U> &m) {
    auto shape = m.shape();
    auto info = m.request();
    T matrix(shape[0], shape[1]);

    vector<U> vec(info.shape[0]*info.shape[1]);

    memcpy(vec.data(), reinterpret_cast<U*>(info.ptr),
        info.size * sizeof(U));

    matrix.load_raw(vec);

    return matrix; 
}

template <class T>
ostream& operator<<(ostream& os, const py::array_t<T> array) {
    auto info = array.request();

    for (int i = 0; i < info.size; i++) {
        os << *static_cast<double*>(info.ptr + i * info.itemsize) << " ";
    }

    return os;
}

template <class T>
ostream& operator<<(ostream& os, const vector<T> vec) {
    for (int i = 0; i < vec.size(); i++) {
        os << vec[i] << " ";
    }

    return os;
}

template <typename T>
vector<T> vfa(py::array_t<T> array) {
    auto buf = array.request();
    T* ptr = static_cast<T*>(buf.ptr);
    return vector<T>(ptr, ptr + buf.size);
}

template <typename T>
py::array_t<T> afv(const vector<T>& vec) {
    return py::array_t<T>(vec.size(), vec.data());
}

PYBIND11_MODULE(lagmodel, m) {
    m.doc() = "Module with functions for the Lag model";
    m.def("find_nearest", &find_nearest<double, vector>,
        "find nearest value in array and report it's index;\n\n"
        "find_nearest(array, value) -> (value, index)", "array"_a, "value"_a);
    m.def("calc_dispfrac", [](py::array_t<double> rad,
                py::array_t<double> rad_area, double rin, double rcor,
                double seedff_norm, double seedff_ind, double heatff_norm,
                double heatff_ind){
            auto [dispfrac, seed_frac_flow, heat_frac_flow] =
                calc_dispfrac<double, vector>(vfa(rad), vfa(rad_area), rin,
                    rcor, seedff_norm, seedff_ind, heatff_norm, heatff_ind);

            return py::make_tuple(afv(dispfrac), afv(seed_frac_flow),
                afv(heat_frac_flow));
        },
        "calculate the disspation fraction;\n\n"
        "calc_dispfrac(rad, rad_area, rin, rcor, seedff_norm, seedff_ind, "
        "heatff_norm, heatff_ind) -> (dispfrac, seef_frac_flow, "
        "heat_frac_flow)", "rad"_a, "rad_area"_a, "rin"_a, "rcor"_a,
        "seedff_norm"_a, "seedff_ind"_a, "heatff_norm"_a, "heatff_ind"_a);
    m.def("calc_illumination_fracs", [](Geometry<vector, double> geomod){
            auto [omega_cor, frad_disktocor, frad_cortodisk] =
                calc_illumination_fracs<double, vector>(geomod);

            return py::make_tuple(afv(omega_cor), afv(frad_disktocor),
                afv(frad_cortodisk));
        },
        "calculate the illumination fractions;\n\n"
        "calc_illumination_fracs(geomod) -> (omega_cor, frad_disktocor, "
        "frad_cortodisk)", "geomod"_a);
    m.def("calc_timing_params", [](py::array_t<double> rad, int i_rsigmax,
                double rcor, int i_rcor, double t_scale,
                tuple<double, double> disk_tau_par,
                tuple<double, double> cor_tau_par, string lor_model,
                py::array_t<double> lor_par, bool dbg){
            auto [tau, lfreq, q_list, rms] =
                calc_timing_params<double, vector>(vfa(rad), i_rsigmax, rcor,
                    i_rcor, t_scale, disk_tau_par, cor_tau_par, lor_model,
                    vfa(lor_par), dbg);

            return py::make_tuple(afv(tau), afv(lfreq), afv(q_list), afv(rms));
        },
        "calculate the timing parameters;\n\n"
        "calc_timing_params(rad, i_rsigmax, rcor, i_rcor, t_scale, "
        "disk_tau_par, cor_tau_par, lor_model, lor_par, dbg=false) -> (tau, "
        "lfreq, q_list, rms)", "rad"_a, "i_rsigmax"_a, "rcor"_a, "i_rcor"_a,
        "t_scale"_a, "disk_tau_par"_a, "cor_tau_par"_a, "lor_model"_a,
        "lor_par"_a, "dbg"_a = false);
    m.def("calc_propagation_params", [](py::array_t<double> rad,
                py::array_t<double> rad_edge, double rcor,
                tuple<double, double> disk_prop_par,
                tuple<double, double> cor_prop_par){
            auto deltau = calc_propagation_params<double, vector>(vfa(rad),
                vfa(rad_edge), rcor, disk_prop_par, cor_prop_par);
            return afv(deltau);
        },
        "calculate the propagation parameters;\n\n"
        "calc_propagation_params(rad, rad_edge, rcor, disk_prop_par, "
        "cor_prop_par) -> (deltau)", "rad"_a, "rad_edge"_a, "rcor"_a,
        "disk_prop_par"_a, "cor_prop_par"_a);
    m.def("calc_radial_time_response", [](py::array_t<double> rad, double rcor,
                py::array_t<double> disp_frac,
                py::array_t<double> disktocor_frac,
                py::array_t<double> cortodisk_frac,
                double thermal_frac, py::array_t<double> seed_frac_flow,
                py::array_t<double> heat_frac_flow){
            auto [ldisk_disp, lseed_disp, lheat, ldisk_rev, lseed_rev] =
                calc_radial_time_response<double, vector>(vfa(rad), rcor,
                    vfa(disp_frac), vfa(disktocor_frac), vfa(cortodisk_frac),
                    thermal_frac, vfa(seed_frac_flow), vfa(heat_frac_flow));
        
            return py::make_tuple(afv(ldisk_disp), afv(lseed_disp), afv(lheat),
                afv(ldisk_rev), afv(lseed_rev));
        },
        "calculate the radial time response;\n\n"
        "calc_radial_time_repsonse(rad, rcor, disp_frac, disktocor_frac, "
        "cortodisk_frac, thermal_frac, seed_frac_flow, heat_frac_flow) -> "
        "(ldisk_disp, lseed_disp, lheat, ldisk_rev, lseed_rev)", "rad"_a,
        "rcor"_a, "disp_frac"_a, "disktocor_frac"_a, "cortodisk_frac"_a,
        "thermal_frac"_a, "seed_frac_flow"_a, "heat_frac_flow"_a);
    m.def("calc_irfs_mono", [](tuple<double, double> gamma_par, double e_seed,
                py::array_t<double> energy, py::array_t<double> ldisk_disp,
                py::array_t<double> lseed_disp, py::array_t<double> lheat,
                py::array_t<double> ldisk_rev, py::array_t<double> lseed_rev,
                bool return_internal){
            auto [gamma_mean, gamma_irf, flux_irf, disk_irf, seed_irf] =
                calc_irfs_mono<double, Array<vector, double>, vector>(gamma_par,
                    e_seed, vfa(energy), vfa(ldisk_disp), vfa(lseed_disp),
                    vfa(lheat), vfa(ldisk_rev), vfa(lseed_rev));

            if (return_internal) {
                return py::make_tuple(gamma_mean, afv(gamma_irf), flux_irf,
                    afv(disk_irf), afv(seed_irf));
            } else {
                return py::make_tuple(gamma_mean, afv(gamma_irf),
                    as_array<double, Array<vector, double>>(flux_irf),
                    afv(disk_irf), afv(seed_irf));
            }
        },
        "caclulate the irfs mono;\n\n"
        "calc_irfs_mono(gama_par, e_seed, energy, ldisk_disp, lseed_disp, "
        "lheat, ldisk_rev, lseed_rev) -> (gamma_mean, gamma_irf, flux_irf, "
        "disk_irf, seed_irf)", "gamma_par"_a, "e_seed"_a, "energy"_a,
        "ldisk_disp"_a, "lseed_disp"_a, "lheat"_a, "ldisk_rev"_a,
        "lseed_rev"_a, "return_internal"_a = false);
    // m.def("calc_irfs_mono", &calc_irfs_mono<double, Nested_Array<list, double>>);
    m.def("calc_disk_band", [](py::array_t<double> disp_frac,
                py::array_t<double> cortodisk_frac,
                py::array_t<double> disktocor_frac,
                py::array_t<double> lseed_disp, py::array_t<double> lheat,
                double thermal_frac, double rcor, py::array_t<double> rad,
                py::array_t<double> rad_area,  tuple<double, double> eband,
                bool kTvar){

            auto [band_frac, kT_rad, ldisk_band, lseed_band, ldiskband_disp, ldiskband_rev] = calc_disk_band<double, vector>(vfa(disp_frac),
                vfa(cortodisk_frac), vfa(disktocor_frac), vfa(lseed_disp),
                vfa(lheat), thermal_frac, rcor, vfa(rad), vfa(rad_area), eband,
                kTvar);

            return py::make_tuple(afv(band_frac), afv(kT_rad), afv(ldisk_band),
                afv(lseed_band), afv(ldiskband_disp), afv(ldiskband_rev));
        },
        "calculate the disk band;\n\n"
        "calc_disk_band(disp_frac, cortodisk_frac, disktocor_frac, lseed_disp, "
        "lheat, thermal_frac, rcor, rad, rad_area, eband, kTvar) -> (band_frac, "
        "kT_rad, ldisk_band, lseed_band, ldiskband_disp, ldiskband_rev)",
        "disp_frac"_a, "cortodisk_frac"_a, "disktocor_frac"_a, "lseed_disp"_a,
        "lheat"_a, "thermal_frac"_a, "rcor"_a, "rad"_a, "rad_area"_a, "eband"_a,
        "kTvar"_a);
    m.def("bb_phflux", &bb_phflux<double>,
        "calculate the black box photon flux;\n\n"
        "bb_phflux(E_min, E_max, kT, nE, kTvar) -> (val)", "E_min"_a, "E_max"_a,
        "kT"_a, "nE"_a, "kTvar"_a);
    m.def("linear_rebin_irf", [](double dt, int i_rsigmax,
                py::array_t<double> irf_nbins,
                py::array_t<double> irf_binedgefrac,
                py::array_t<double> input_irf, py::array_t<double> deltau_scale,
                int nirf){
            auto [rebinned_irf, flux_outer] =
                linear_rebin_irf<double, vector>(dt, i_rsigmax, vfa(irf_nbins),
                    vfa(irf_binedgefrac), vfa(input_irf), vfa(deltau_scale),
                    nirf);

            return py::make_tuple(afv(rebinned_irf), flux_outer);
        },
        "rebin to linear scale;\n\n"
        "linear_rebin_irf(dt, i_rsigmax, irf_nbins, irf_binedgefrac, "
        "input_irf, deltau_scale, nirf) -> (rebinned_irf, flux_outer)", "dt"_a,
        "i_rsigmax"_a, "irf_nbins"_a, "irf_binedgefrac"_a, "input_irf"_a,
        "deltau_scale"_a, "nirf"_a);
    m.def("calc_cross_psd", [](py::array_t<double> freq, double dt,
                py::array_t<double> ci_irf, py::array_t<double> ref_irf,
                py::array_t<double> irf_nbins, py::array_t<double> deltau_scale,
                int i_rsigmax, py::array_t<double> f_pk, py::array_t<double> q,
                py::array_t<double> rms){
            auto [wt_cross_spec, wt_pow_spec_ci, wt_pow_spec_ref, mod_sig_psd] =
                calc_cross_psd<double, vector>(vfa(freq), dt, vfa(ci_irf),
                    vfa(ref_irf), vfa(irf_nbins), vfa(deltau_scale), i_rsigmax,
                    vfa(f_pk), vfa(q), vfa(rms));
        
            return py::make_tuple(afv(wt_cross_spec), afv(wt_pow_spec_ci),
                afv(wt_pow_spec_ref), afv(mod_sig_psd));
        },
        "calculate the cross power spectra\n\n"
        "calc_cross_psd(freq, dt, ci_irf, ref_irf, irf_nbins, deltau_scale, "
        "i_rsigmax, f_pk, q, rms) -> (wt_cross_spec, wt_pow_spec_ci, "
        "wt_pow_spec_ref, mod_sig_psd)", "freq"_a, "dt"_a, "ci_irf"_a,
        "ref_irf"_a, "irf_nbins"_a, "deltau_scale"_a, "i_rsigmax"_a, "f_pk"_a,
        "q"_a, "rms"_a);
    m.def("lorentz_q", [](py::array_t<double> f, py::array_t<double> f_pk,
                py::array_t<double> q, py::array_t<double> rms,
                bool return_internal){
            auto lorentz =
                lorentz_q<double, Array<vector, double>, vector>(vfa(f),
                    vfa(f_pk), vfa(q), vfa(rms));

            variant<py::array_t<double>, Array<vector, double>> result;
            if (return_internal) {
                result = lorentz;
            } else {
                result = as_array<double, Array<vector, double>>(lorentz);
            }

            return result;
        },
        "f"_a, "f_pk"_a, "q"_a , "rms"_a, "return_internal"_a = false);
    m.def("lorentz_q", [](py::array_t<double> f, double f_pk, double q,
                double rms){
            auto lorentz = lorentz_q<double, vector>(vfa(f), f_pk, q, rms);

            return afv(lorentz);
        },
        "calculates the lorentzian function;\n\n"
        "lorentz_q(f, f_pk, q, rms) -> (lorentz)", "f"_a, "f_pk"_a, "q"_a,
        "rms"_a);
    m.def("calculate_stprod_mono", [](int nirf_mult, py::array_t<double> energy,
                py::array_t<int> encomb, py::array_t<double> flux_irf,
                py::array_t<double> disk_irf, py::array_t<double> gamma_irf,
                py::array_t<double> deltau, double min_deltau_frac,
                int i_rsigmax, py::array_t<double> lfreq, py::array_t<double> q,
                py::array_t<double> rms, double t_scale, bool dbg) {
            assert(encomb.ndim() == 2 && flux_irf.ndim() == 2);

            auto en = as_matrix<Array<vector, int>, int>(encomb);
            auto flux = as_matrix<Array<vector, double>, double>(flux_irf);

            auto [freq, phlag, tlag, psd_ci, psd_ref, mod_sig_psd, irf_nbins,
                irf_binedgefrac, deltau_scale, dt, nirf, ci_irf, ci_mean,
                ci_outer] = calculate_stprod_mono<double, Array<vector, double>,
                    Array<vector, int>>(nirf_mult, vfa(energy), en, flux,
                        vfa(disk_irf), vfa(gamma_irf), vfa(deltau),
                        min_deltau_frac, i_rsigmax, vfa(lfreq), vfa(q), vfa(rms),
                        t_scale, dbg);
            
            return py::make_tuple(afv(freq),
                as_array<double, Array<vector, double>>(phlag),
                as_array<double, Array<vector, double>>(tlag),
                as_array<double, Array<vector, double>>(psd_ci),
                as_array<double, Array<vector, double>>(psd_ref),
                afv(mod_sig_psd), afv(irf_nbins), afv(irf_binedgefrac),
                afv(deltau_scale), dt, nirf,
                as_array<double, Array<vector, double>>(ci_irf), afv(ci_mean),
                afv(ci_outer));
        },
        "calculates mono-energetic spectral-timing products;\n\n"
        "calculate_stprod_mono(nirf_mult, energy, encomb, flux_irf, disk_irf, "
        "gamma_irf, deltau, min_deltau_frac, i_rsigmax, lfreq, q, rms, "
        "t_scale, dbg=false) -> (freq, phlag, tlag, psd_ci, psd_ref, "
        "mod_sig_psd,irf_nbins, irf_binedgefrac, deltau_scale, dt, nirf, "
        "ci_irf, ci_mean, ci_outer)", "nirf_mult"_a, "energy"_a, "encomb"_a,
        "flux_irf"_a, "disk_irf"_a, "gamma_irf"_a, "deltau"_a,
        "min_deltau_frac"_a, "i_rsigmax"_a, "lfreq"_a, "q"_a, "rms"_a,
        "t_scale"_a, "dbg"_a = false);

    auto internal = m.def_submodule("internal_class");

    typedef Matrix<vector, double> matrix_0;

    py::class_<matrix_0>(internal, "Matrix", py::buffer_protocol())
        .def(py::init<>())
        .def("get_size", &matrix_0::get_size);

    typedef Nested_Array<vector, double> matrix_1;

    py::class_<matrix_1, matrix_0>(internal, "Nested_Array",
            py::buffer_protocol())
        .def(py::init<list<double>, list<double>>())
        .def(py::init<vector<double>, vector<double>>())
        .def(py::init<matrix_1>())
        .def(py::init<int, int>())
        .def("get_element", py::overload_cast<int, int>(&matrix_1::get_element))
        .def("get_element", py::overload_cast<int, int>(&matrix_1::get_element,
            py::const_))
        .def("set_element", &matrix_1::set_element)
        .def("get_matrix", &matrix_1::get_matrix)
        .def_buffer([](matrix_1 &m) -> py::buffer_info {
            auto [row, col] = m.get_size();
            auto ptr_val = m.get_matrix().data();

            return py::buffer_info(
                &ptr_val, sizeof(double),
                py::format_descriptor<double>::format(), 2, {row, col},
                {sizeof(double) * col, sizeof(double)}
            );
        });

    typedef Array<vector, double> matrix_2;

    py::class_<matrix_2, matrix_0>(internal, "Array", py::buffer_protocol())
        .def(py::init<list<double>, list<double>>())
        .def(py::init<vector<double>, vector<double>>())
        .def(py::init<matrix_2>())
        .def(py::init<int, int>())
        .def(py::init([](py::buffer b) {
            py::buffer_info info = b.request();
            if (info.format != py::format_descriptor<double>::format()) {
                // cout << info.format << endl;
                // cout << py::format_descriptor<double>::format() << endl;
                throw std::runtime_error("Incompatible format: expected a "
                "double array!");
            }
            if (info.ndim != 2) {
                throw std::runtime_error("Incompatible buffer dimension!");
            }

            matrix_2 m(info.shape[0], info.shape[1]);

            vector<double> vec(info.shape[0]*info.shape[1]);

            memcpy(vec.data(), reinterpret_cast<double *>(info.ptr),
                info.size * sizeof(double));

            m.load_raw(vec);

            return m;
        }))
        .def("get_element", py::overload_cast<int, int>(&matrix_2::get_element))
        .def("get_element", py::overload_cast<int, int>(&matrix_2::get_element,
            py::const_))
        .def("set_element", &matrix_2::set_element)
        .def("get_matrix", &matrix_2::get_matrix)
        .def_buffer([](matrix_2 &m) -> py::buffer_info {
            auto [row, col] = m.get_size();
            auto ptr_val = m.get_matrix().data();

            return py::buffer_info(
                ptr_val, sizeof(double),
                py::format_descriptor<double>::format(), 2, {row, col},
                {sizeof(double) * col, sizeof(double)}
            );
        });

    internal.def("test", [](matrix_2 m) {
        auto [row, col] = m.get_size();

        for (int i = 0; i < row; i++) {
            for (int j = 0; j < col; j++) {
                cout << m.get_element(i, j) << " ";
            } cout << endl;
        } 
    });
    internal.def("test", [internal](py::array_t<double> m) {
        assert(m.ndim() == 2);
        auto shape = m.shape();

        for (int i = 0; i < shape[0]; i++) {
            for (int j = 0; j < shape[1]; j++) {
                cout << m.at(i, j) << " ";
            } cout << endl;
        }

        cout << endl << "----" << endl;
        as_matrix<Array<vector, double>, double>(m);
        // auto a = internal.attr("Array")(m.request());
        // auto a = Array<vector, double>(m.shape()[0], m.shape()[1]);
        // a.load_raw(m.data());
    });

    auto geo = m.def_submodule("geometries", "Geometry classes for the lagmodel"
        " module");

    typedef Geometry<vector, double> geo_base;

    py::class_<geo_base>(geo, "Geometry")
        .def(py::init<vector<double>, vector<double>, double>())
        .def("get_rcor", &geo_base::get_rcor)
        .def("get_omega_cor", &geo_base::get_omega_cor)
        .def("get_frad_disktocor", &geo_base::get_frad_disktocor)
        .def("get_frad_cortodisk", &geo_base::get_frad_cortodisk)
        .def("get_omega_cor_area", &geo_base::get_omega_cor_area)
        .def("get_frad_disktocor_area", &geo_base::get_frad_disktocor_area)
        .def("get_frad_cortodisk_area", &geo_base::get_frad_cortodisk_area);

    typedef AN_Sphere<vector, double> geo_ansphere;

    py::class_<geo_ansphere, geo_base>(geo, "AN_Sphere")
        .def(py::init<vector<double>, vector<double>, double>());

    typedef BKNpow_Emiss<vector, double> geo_bknpowemiss;

    py::class_<geo_bknpowemiss, geo_base>(geo, "BKNpow_Emiss")
        .def(py::init<vector<double>, vector<double>, double, double, double,
            double, double, double, double>());

    typedef Cylinder<vector, double> geo_cylinder;

    py::class_<geo_cylinder, geo_base>(geo, "Cylinder")
        .def(py::init<vector<double>, vector<double>, double, double, int,
            int>())
        .def(py::init<vector<double>, vector<double>, double, double, double,
            double>());

    typedef Inv_Cone<vector, double> geo_invcone;

    py::class_<geo_invcone, geo_base>(geo, "Inv_Cone")
        .def(py::init<vector<double>, vector<double>, double, double, double,
            int, int>())
        .def(py::init<vector<double>, vector<double>, double, double, double,
            double, double>());

    typedef Sphere<vector, double> geo_sphere;

    py::class_<geo_sphere, geo_base>(geo, "Sphere")
        .def(py::init<vector<double>, vector<double>, double, int, int>())
        .def(py::init<vector<double>, vector<double>, double, double, double>());
}