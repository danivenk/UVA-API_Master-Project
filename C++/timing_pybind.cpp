
// include libraries
#include <variant>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11_json/pybind11_json.hpp>
#include "timing_functions.h"

// define some name spaces
namespace py = pybind11;
using namespace py::literals;

/**
 * @brief Convert a matrix class to a numpy array.
 * 
 * @tparam T is the value type
 * @tparam U is the matrix type
 * @param m is the matrix to be converted
 * @return returns the matrix as numpy array
 */
template<class T, class U>
py::array_t<T> as_array(U &m) {

    // get the number of rows and columns, the pointer to the data
    auto [row, col] = m.get_size();
    auto ptr_val = m.get_matrix().data();

    // make a copy of the data
    auto data_copy = std::make_unique<T[]>(row * col);
    memcpy(data_copy.get(), ptr_val, row * col * sizeof(T));

    // return it as a numpy array
    return py::array_t<T>({row, col}, {sizeof(T) * col, sizeof(T)},
        data_copy.get());
}

/**
 * @brief Convert a numpy array matrix to a matrix class.
 * 
 * @tparam T is the matrix type
 * @tparam U is the value type
 * @param m is the numpy matrix to be converted
 * @return returns the matrix
 */
template<class T, class U>
T as_matrix(py::array_t<U> &m) {

    // get te shape and info and make an empty matrix
    auto shape = m.shape();
    auto info = m.request();
    T matrix(shape[0], shape[1]);

    // make a vector of the numpy data and copy to it
    vector<U> vec(info.shape[0]*info.shape[1]);

    memcpy(vec.data(), reinterpret_cast<U*>(info.ptr),
        info.size * sizeof(U));

    // load this vector into the matrix
    matrix.load_raw(vec);

    // return the matrix
    return matrix; 
}

/**
 * @brief Stream operator for a numpy array.
 * 
 * @tparam T is the value type
 * @param os is the stream to output to
 * @param array is the array to be streamed
 * @return returns the output stream
 */
template <class T>
ostream& operator<<(ostream& os, const py::array_t<T> array) {

    // get the data
    auto info = array.request();

    // print it all to the output stream
    for (int i = 0; i < info.size; i++) {
        os << *static_cast<double*>(info.ptr + i * info.itemsize) << " ";
    }

    // return the output stream
    return os;
}

/**
 * @brief Stream operator for a vector
 * 
 * @tparam T is the value type
 * @param os is the stream to outpu to
 * @param vec is the vector to be streamed
 * @return returns the output stream
 */
template <class T>
ostream& operator<<(ostream& os, const vector<T> vec) {

    // loop over the vector elements and output to the stream
    for (int i = 0; i < vec.size(); i++) {
        os << vec[i] << " ";
    }

    // return the output stream
    return os;
}

/**
 * @brief Create a vector from a numpy array.
 * 
 * @tparam T is the value type 
 * @param array is the numpy array to be converted
 * @return returns the converted numpy array as a vector
 */
template <typename T>
vector<T> vfa(py::array_t<T> array) {

    // get the data and cast the pointer to the data to something usable
    auto buf = array.request();
    T* ptr = static_cast<T*>(buf.ptr);

    // returns the converted numpy array as a vector
    return vector<T>(ptr, ptr + buf.size);
}

/**
 * @brief Create a numpy array from a vector
 * 
 * @tparam T is the value type
 * @param vec is the vector to be converted
 * @return returns the vector as a converted numpy array
 */
template <typename T>
py::array_t<T> afv(const vector<T>& vec) {

    // returns the vector as a numpy array
    return py::array_t<T>(vec.size(), vec.data());
}

// create the module
PYBIND11_MODULE(timing_lagmodel, m) {

    // describe the module and bind each function and classes
    m.doc() = "Module with functions for the Lag model";
    init_json();
    m.def("find_nearest", &find_nearest<double, vector>,
        "Find the nearest value in an array.\n\n:param Array[float] array: is "
        "the array to be searched\n:param float value: the value to find the "
        "nearest value for\n:return: returns the closest value and it's index\n"
        ":rtype: Tuple[float, int]",
        "array"_a, "value"_a);
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
        "Calculate the dissipation for each radial bin normalized by the total "
        "dissipation for all radial bins. It also calculates the heat and seed "
        "photon fraction of the dissipated power.\n\nTo allow for unconstrained"
        " heating/cooling the `seedff_norm` can be set as bigger than 0. If it "
        "is set to less than 0 to force the total fraction to be unity.\nWhen "
        "`seedff_norm` is set to 0 no seed luminosity is produced inside the "
        "corona.\n\n:param Array[float] rad: radial bins\n:param Array[float] "
        "rad_area: radial area bins\n:param float rin: a radius inside the "
        "corona\n:param float rcor: radius of the corona\n:param float "
        "seedff_norm: normalization factor for the seed fraction flow\n:param "
        "float seedff_ind: exponent for the seed fraction flow\n:param float "
        "heatff_norm: normalization factor for the heat fraction flow\n:param "
        "float heatff_ind: exponent for the heat fraction flow\n:return: "
        "returns the dissipation fraction, seed fraction flow, heat fraction "
        "flow\n:rtype: Tuple[Array[float], Array[float], Array[float]]",
        "rad"_a, "rad_area"_a, "rin"_a, "rcor"_a, "seedff_norm"_a,
        "seedff_ind"_a, "heatff_norm"_a, "heatff_ind"_a);
    m.def("calc_illumination_fracs", [](Geometry<vector, double> geomod){
            auto [omega_cor, frad_disktocor, frad_cortodisk] =
                calc_illumination_fracs<double, vector>(geomod);

            return py::make_tuple(afv(omega_cor), afv(frad_disktocor),
                afv(frad_cortodisk));
        },
        "Calculates the solid angle & illumination fractions. The illumination "
        "fractions are the fraction of emission that goes from the corona to "
        "the disk and the emission that goes from the disk to the corona.\n\n"
        ":param Geometry geomod: is the geometry to calculate the illumination "
        "fractions and solid angle for\n:return: returns the solid angle and "
        "illumination fractions for the given geometry\n:rtype: "
        "Tuple[Array[float], Array[float], Array[float]]",
        "geomod"_a);
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
        "Calculates the timing parameters for a given Lorentzian. It these "
        "include the assigned signal power spectrum for each bin from the "
        "Lorentzian, the Lorentzian frequency and the matching Lorentzian RMS "
        "value.\n\nFor this calculation 4 different Lorentzian models exits:\n-"
        " continuous (Lorentzian signals at all radial bins) ~\nLorentzian "
        "model parameters: [q, total_rms_in_disk, total_rms_in_corona];\n- "
        "multi_frequency (multiple Lorentzian parameters assigned dependent on "
        "tau) ~\nLorentzian model parameters: [peak_frequency_1, q_1, rms_1, "
        "peak_frequency_2, q_2, rms_2, ...] where len(lor_par) % 3 == 0;\n- "
        "multi_radius (all the parameters set dependent on the radial bins) ~\n"
        "Lorentzian model parameters: [radius_1, q_1, rms_1, radius_2, q_2, "
        "rms_2, ...] where len(lor_par) % 3 == 0;\n- multi_radius_frequency (is"
        " a combination of the two previous) ~\nLorentzian model parameters: "
        "[radius_1, peak_frequency_1, q_1, rms_1, radius_2, peak_frequency_2, "
        "q_2, rms_2, ...] where len(lor_par) % 4 == 0;\n\n:param Array[float] "
        "rad: radial bins\n:param int i_rsigmax: index of the nearest value for"
        " the radial bin for signal max\n:param float rcor: radius of the "
        "corona\n:param int i_rcor: index of the nearest radial bin to the "
        "coronal radius\n:param t_scale: scaling factor to convert the "
        "propagation delay to seconds\n:param Tuple[float, float] disk_tau_par:"
        " parameters for the variability time-scale in the disk\n:param "
        "Tuple[float, float] cor_tau_par: parameters for the variability "
        "time-scale in the corona\n:param str lor_model: Lorenzian model to use"
        "\n:param Array[float] lor_par: parameters for the chosen Lorentzian\n"
        ":param bool dbg: print out signal radius (default: False)\n:return: "
        "returns the variable time scale, Lorentzian frequency, quality factor,"
        " rms values\n:rtype: Tuple[Array[float], Array[float], Array[float], "
        "Array[float]]",
        "rad"_a, "i_rsigmax"_a, "rcor"_a, "i_rcor"_a, "t_scale"_a,
        "disk_tau_par"_a, "cor_tau_par"_a, "lor_model"_a, "lor_par"_a,
        "dbg"_a = false);
    m.def("calc_propagation_params", [](py::array_t<double> rad,
                py::array_t<double> rad_edge, double rcor,
                tuple<double, double> disk_prop_par,
                tuple<double, double> cor_prop_par){
            auto deltau = calc_propagation_params<double, vector>(vfa(rad),
                vfa(rad_edge), rcor, disk_prop_par, cor_prop_par);
            return afv(deltau);
        },
        "Calculates the propagation delay across the radial bins according to "
        "the given propagation parameters.\n\n:param Array[float] rad: radial "
        "bins\n:param Array[float] rad_edge: edges of the radial bins\n:param "
        "float rcor: radius of the corna\n:param Tuple[float, float] "
        "disk_prop_par: propagation parameters for the disk\n:param "
        "Tuple[float, float] cor_prop_par: propagation parameters for the "
        "corona\n:return: returns the propagation delay\n:rtype: Array[float]",
        "rad"_a, "rad_edge"_a, "rcor"_a, "disk_prop_par"_a, "cor_prop_par"_a);
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
        "Calculate the disk, seed & heat impulse responses due to the "
        "disspation from each bin.\n\n:param Array[float] rad: radial bins\n"
        ":param float rcor: radius of the corona\n:param Array[float] "
        "disp_frac: dissipation fraction\n:param disktocor_frac: luminosity "
        "fraction from the disk to the corona\n:param Array[float] "
        "cortodisk_frac: luminosity fraction from the corona to the disk\n"
        ":param float thermal_frac: fraction of luminosity not turned into "
        "thermal energy\n:param Array[float] seed_frac_flow: fraction of seed "
        "photons in the flow\n:param Array[float] heat_frac_flow: fraction of "
        "heat photons in the flow\n:return: returns the luminosities for the "
        "disk disspation, seed disspation, heating luminosity, disk "
        "reverberation, seed reverberation\n:rtype: Tuple[Array[float], "
        "Array[float], Array[float], Array[float]]",
        "rad"_a, "rcor"_a, "disp_frac"_a, "disktocor_frac"_a,
        "cortodisk_frac"_a, "thermal_frac"_a, "seed_frac_flow"_a,
        "heat_frac_flow"_a);
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
        "Calculates the impulse response functions (IRFs) for the monochromatic"
        " photon flux, photon index, total disk emission and total seed "
        "emission. For the total emissions they are summed over both the "
        "dissipation and reverberation emissions. The IRFs are calculated "
        "depending on the given luminosities.\n\n:param Tuple[float, float] "
        "gamma_par: parameters for the photon index\n:param float e_seed: "
        "monochromatic approximation for the seed energy\n:param Array[float] "
        "energy: energy value to use\n:param Array[float] ldisk_disp: "
        "luminosity for the disk dissipation\n:param Array[float] lseed_disp: "
        "luminosity for the seed dissipation\n:param Array[float] lheat: "
        "heating luminosity\n:param Array[float] ldisk_rev: luminosity for the "
        "disk reverberation\n:param Array[float] lseed_rev: luminosity for the "
        "seed reverberation\n:param return_internal: false = returns Array, "
        "true = returns internal Matrix (default: False)\n:return: returns the "
        "photon index mean, impulse response function for the photon index, photon flux, total disk emission, total seed emission\n:rtype: Tuple[float, Array[float], Matrix[float], Array[float], Array[float]]",
        "gamma_par"_a, "e_seed"_a, "energy"_a, "ldisk_disp"_a, "lseed_disp"_a,
        "lheat"_a, "ldisk_rev"_a, "lseed_rev"_a, "return_internal"_a = false);
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
        "Calculate the luminosities emitted for a given black body energy band "
        "as well as the fraction of photons.\n\n:param Array[float] disp_frac: "
        "dissipation fraction\n:param Array[float] cortodisk_frac: luminosity "
        "fraction from the disk to the corona\n:param Array[float] "
        "disktocor_frac: luminosity fraction from the corona to the disk\n"
        ":param Array[float] lseed_disp: luminosity of the seed dissipation\n"
        ":param Array[float] lheat: heating luminosity\n:param float "
        "thermal_frac: fraction of luminosity not turned into thermal energy\n"
        ":param float rcor: radius of the corona\n:param Array[float] rad: "
        "radial bins\n:param Array[float] rad_area: radial area bins\n:param "
        "Tuple[float, float] eband: the energy band to use\n:param bool kTvar: "
        "true = non-variable black body; false = variable black body\n:return: "
        "returns the band fraction, temperature radial bins, luminosities for "
        "this energy band: disk, seed, disk dissipation, disk reverberation\n"
        ":rtype: Tuple[Array[float], Array[float], Array[float], Array[float], "
        "Array[float], Array[float]]",
        "disp_frac"_a, "cortodisk_frac"_a, "disktocor_frac"_a, "lseed_disp"_a,
        "lheat"_a, "thermal_frac"_a, "rcor"_a, "rad"_a, "rad_area"_a, "eband"_a,
        "kTvar"_a);
    m.def("bb_phflux", &bb_phflux<double>,
        "Calculates the black body photon flux for the given values.\n\n:param "
        "float E_min: is the minimum energy value to probe at\n:param float "
        "E_max: is the maximum energy value to probe at\n:param float kT: is "
        "the temperature\n:param int nE: is the number of energy bins\n:param "
        "bool kTvar: true = non-variable black body; false = variable black "
        "body\n:return: returns the black body photon flux\n:rtype: float",
        "E_min"_a, "E_max"_a, "kT"_a, "nE"_a, "kTvar"_a);
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
        "Linearly rebin the impulse response function (IRF) given. The output "
        "IRF are only filled to the maximum signal radius to save memory.\n\n"
        ":param float dt: is the time step\n:param int i_rsigmax: index of the "
        "nearest value for the radial bin for signal max\n:param Array[float] "
        "irf_nbins: number of bins for the impulse response function\n:param "
        "Array[float] irf_binedgefrac: fraction at the edge of the bin for the "
        "IRF\n:param Array[float] input_irf: input impulse response function\n"
        ":param Array[float] deltau_scale: scale for the propagation delay\n"
        ":param int nirf: number of impulse response functions\n:return: "
        "returns the rebinned impulse response function & the flux component "
        "outside the maximum singal radius\n:rtype: Tuple[Array[float], float]",
        "dt"_a, "i_rsigmax"_a, "irf_nbins"_a, "irf_binedgefrac"_a,
        "input_irf"_a, "deltau_scale"_a, "nirf"_a);
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
        "Calculate the cross and power spectrum for the given impulse response function (IRF) and find driving signal power spectrum density (PSD)\n\n"
        ":param Array[float] freq: frequencies\n:param float dt: time step\n"
        ":param Array[float] ci_irf: channel of interest impulse response "
        "function\n:param Array[float] ref_irf: reference impulse response "
        "function\n:param Array[float] irf_nbins: number of bins for the "
        "impulse response function\n:param Array[float] deltau_scale: scale for"
        " the propagation delay\n:param int i_rsigmax: index of the nearest "
        "value for the radial bin for signal max\n:param Array[float] f_pk: "
        "peak frequencies\n:param Array[float] q: quality factor\n:param "
        "Array[float] rms: root mean square\n:return: returns the weighted "
        "cross spectrum, signal power spectrum, reference power spectrum and "
        "the signal power spectrum density\n:rtype: Tuple[Array[complex[float]]"
        ", Array[float], Array[float], Array[float]]",
        "freq"_a, "dt"_a, "ci_irf"_a, "ref_irf"_a, "irf_nbins"_a,
        "deltau_scale"_a, "i_rsigmax"_a, "f_pk"_a, "q"_a, "rms"_a);
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
        "Calculate the Lorentzian Signal power density.\n\n:param Array[float] "
        "f: frequencies\n:param Array[float] f_pk: peak frequencies\n:param "
        "Array[float] q: quality value\n:param Array[float] rms: root mean "
        "square\n:param return_internal: false = returns Array, true = returns "
        "internal Matrix (default: False)\n:return: returns matrix of "
        "Lorentzian Signal power densities\n:rtype: Array[float] or "
        "Matrix[float]", 
        "f"_a, "f_pk"_a, "q"_a , "rms"_a, "return_internal"_a = false);
    m.def("lorentz_q", [](py::array_t<double> f, double f_pk, double q,
                double rms){
            auto lorentz = lorentz_q<double, vector>(vfa(f), f_pk, q, rms);

            return afv(lorentz);
        },
        "Calculate the Lorentzian Signal power density\n\n:param Array[float] "
        "f: frequencies\n:param float f_pk: peak frequency\n:param float q: "
        "quality value\n:param float rms: root mean square\n:return: returns "
        "array of Lorentzian Signal power densities\n:rtype: Array[float]",
        "f"_a, "f_pk"_a, "q"_a, "rms"_a);
    m.def("calculate_stprod_mono", [](int nirf_mult, py::array_t<double> energy,
                py::array_t<int> encomb, py::array_t<double> flux_irf,
                py::array_t<double> disk_irf, py::array_t<double> gamma_irf,
                py::array_t<double> deltau, double min_deltau_frac,
                int i_rsigmax, py::array_t<double> lfreq, py::array_t<double> q,
                py::array_t<double> rms, double t_scale, bool dbg,
                bool return_internal) {
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
            
            if (return_internal) {
                return py::make_tuple(afv(freq), phlag, tlag, psd_ci, psd_ref,
                    afv(mod_sig_psd), afv(irf_nbins), afv(irf_binedgefrac),
                    afv(deltau_scale), dt, nirf, ci_irf, afv(ci_mean),
                    afv(ci_outer));
            } else {
                return py::make_tuple(afv(freq), as_array<double, Array<vector,
                    double>>(phlag), as_array<double, Array<vector,
                    double>>(tlag), as_array<double, Array<vector,
                    double>>(psd_ci), as_array<double, Array<vector,
                    double>>(psd_ref), afv(mod_sig_psd), afv(irf_nbins),
                    afv(irf_binedgefrac), afv(deltau_scale), dt, nirf,
                    as_array<double, Array<vector, double>>(ci_irf),
                    afv(ci_mean), afv(ci_outer));
            }
        },
        "Calculates the lag frequency, the power spectra for the disk, seed and"
        " the power-law emission from the different bands.\n\n:param int "
        "nirf_mult: is a multiplication factor for the IRF bins\n:param "
        "Array[float] energy: energy values to use\n:param Matrix[int] encomb: "
        "indices for the channel of interest and reference band\n:param "
        "Matrix[float] flux_irf: impulse response function for the flux\n:param"
        " Array[float] disk_irf: impulse response function for the disk\n:param"
        " Array[float] gamma_irf: impulse response function for the photon "
        "index\n:param Array[float] deltau: propagation delay\n:param float "
        "min_deltau_frac: minimum fraction for the propagation delay\n:param "
        "int i_rsigmax: index for the maximum signal radius\n:param "
        "Array[float] lfreq: Lorentzian frequencies\n:param Array[float] q: "
        "quality values\n:param Array[float] rms: the root mean square values\n"
        ":param float t_scale: scaling factor to convert the propagation delay "
        "to seconds\n:param bool dbg: print or not print out some values to the"
        " console (default: False)\n:param return_internal: false = returns "
        "Array, true = returns internal Matrix (default: False)\n:return: "
        "returns the frequencies, photon and time lags, power spectrum densities for the channel of interest and reference IRFs, the signal power spectrum "
        "density, number of bins for the impulse response function, fraction at"
        " the edge of the bin for the impulse response function, the scaled "
        "propagation delays, the time steps, the number of impulse response "
        "function bins, the channel of interest impulse response function, mean"
        " and outer values\n:rtype: Tuple[Array[float], Matrix[float], "
        "Matrix[float], Matrix[float], Matrix[float], Array[float], "
        "Array[float], Array[float], Array[float], Array[float], float, int, "
        "Matrix[float], Array[float], Array[float]]",
        "nirf_mult"_a, "energy"_a, "encomb"_a, "flux_irf"_a, "disk_irf"_a,
        "gamma_irf"_a, "deltau"_a, "min_deltau_frac"_a, "i_rsigmax"_a,
        "lfreq"_a, "q"_a, "rms"_a, "t_scale"_a, "dbg"_a = false,
        "return_internal"_a = false);
    m.def("get_times",  []{ return times_Cpp; },
        "Get the times JSON.\n\n:return: returns the times JSON\n"
        ":rtype: Dict[str, List[int]]");
    m.def("add_times", [](string key, float val){
            times_Cpp[key].push_back(val);
        },"Add to the times JSON.\n\n:param str key: Key to add at\n:param int "
        "val: Value to add", "key"_a, "val"_a);
    m.def("reset_times", []{ times_Cpp.clear(); }, "Clears the times JSON");

    // create the internal class submodule
    auto internal = m.def_submodule("internal_class");

    // define the base class type
    typedef Matrix<vector, double> matrix_0;

    // define the base class and it's methods
    py::class_<matrix_0>(internal, "Matrix",
            "Defines a base class Matrix with a size.",
            py::buffer_protocol())
        .def(py::init<>(),
            "Inizializes a matrix.")
        .def("get_size", &matrix_0::get_size,
            "Get the size of the matrix\n\n:return: returns the size of the "
            "matrix\n:rtype: Tuple[int, int]");

    // define the first child class type
    typedef Nested_Array<vector, double> matrix_1;

    // define the first child class and it's methods
    py::class_<matrix_1, matrix_0>(internal, "Nested_Array",
            "Defines a matrix saved as a nested array and is based on the "
            "Matrix class", 
            py::buffer_protocol())
        .def(py::init<list<double>, list<double>>(),
            "Initializes a matrix from two lists that are then multiplied with "
            "each other elementwise.\n\nExample:\nNested_Array([a, b, c], [d, "
            "e, f]) => [[ad, ae, af], [bd, be, bf], [cd, ce, cf]]\n\n:param "
            "List[float] arr1: array 1\n:param List[float] arr2: array 2",
            "arr1"_a, "arr2"_a)
        .def(py::init<vector<double>, vector<double>>(),
            "Initializes a matrix from two lists that are then multiplied with "
            "each other elementwise.\n\nExample:\nNested_Array([a, b, c], [d, "
            "e, f]) => [[ad, ae, af], [bd, be, bf], [cd, ce, cf]]\n\n:param "
            "List[float] arr1: array 1\n:param List[float] arr2: array 2",
            "arr1"_a, "arr2"_a)
        .def(py::init<matrix_1>(),
            "Initializes a matrix from a different matrix copying it's "
            "contents.\n\n:param Nested_Array other: the matrix to copy from",
            "other"_a)
        .def(py::init<int, int>(),
            "Initializes a matrix from 2 sizes.\n\n:param int size1: first size"
            "\n:param int size2: second size",
            "size1"_a, "size2"_a)
        .def("get_element", py::overload_cast<int, int>(&matrix_1::get_element),
            "Gets an element at a specific (i,j) location in the matrix.\n\n"
            ":param int i: location for the first index\n:param int j: location"
            " for the second index\n:return: returns the element at location "
            "(i, j) in the matrix\n:rtype: float",
            "i"_a, "j"_a)
        .def("get_element", py::overload_cast<int, int>(&matrix_1::get_element,
            py::const_),
            "Gets an element at a specific (i,j) location in the matrix.\n\n"
            ":param int i: location for the first index\n:param int j: location"
            " for the second index\n:return: returns the element at location "
            "(i, j) in the matrix\n:rtype: float (const)",
            "i"_a, "j"_a)
        .def("set_element", &matrix_1::set_element,
            "Sets an element at a specifi (i, j) location in the matrix.\n\n"
            ":param int i: location for the first index\n:param int j: location"
            " for the second index")
        .def("get_matrix", &matrix_1::get_matrix,
            "Get the internal matrix.\n\n:return: returns the internally stored"
            " matrix\n:rtype: List[List[float]]")
        .def_buffer([](matrix_1 &m) -> py::buffer_info {

            // get sizes and data
            auto [row, col] = m.get_size();
            auto ptr_val = m.get_matrix().data();

            // return buffer from this info
            return py::buffer_info(
                &ptr_val, sizeof(double),
                py::format_descriptor<double>::format(), 2, {row, col},
                {sizeof(double) * col, sizeof(double)}
            );
        });

    // define the second child class type
    typedef Array<vector, double> matrix_2;

    // define the second child class and it's methods
    py::class_<matrix_2, matrix_0>(internal, "Array",
            "Defines a matrix saved as a 1d array and is based on the "
            "Matrix class",
            py::buffer_protocol())
        .def(py::init<list<double>, list<double>>(),
            "Initializes a matrix from two lists that are then multiplied with "
            "each other elementwise.\n\nExample:\nNested_Array([a, b, c], [d, "
            "e, f]) => [[ad, ae, af], [bd, be, bf], [cd, ce, cf]]\n\n:param "
            "List[float] arr1: array 1\n:param List[float] arr2: array 2",
            "arr1"_a, "arr2"_a)
        .def(py::init<vector<double>, vector<double>>(),
            "Initializes a matrix from two lists that are then multiplied with "
            "each other elementwise.\n\nExample:\nNested_Array([a, b, c], [d, "
            "e, f]) => [[ad, ae, af], [bd, be, bf], [cd, ce, cf]]\n\n:param "
            "List[float] arr1: array 1\n:param List[float] arr2: array 2",
            "arr1"_a, "arr2"_a)
        .def(py::init<matrix_2>(),
            "Initializes a matrix from a different matrix copying it's "
            "contents.\n\n:param Nested_Array other: the matrix to copy from",
            "other"_a)
        .def(py::init<int, int>(),
            "Initializes a matrix from 2 sizes.\n\n:param int size1: first size"
            "\n:param int size2: second size",
            "size1"_a, "size2"_a)
        .def(py::init([](py::buffer b) {

            // get buffer info
            py::buffer_info info = b.request();

            // make sure the dimensions and format is correct
            if (info.format != py::format_descriptor<double>::format()) {
                throw std::runtime_error("Incompatible format: expected a "
                "double array!");
            }
            if (info.ndim != 2) {
                throw std::runtime_error("Incompatible buffer dimension!");
            }

            // create a matrix in the correct shape
            matrix_2 m(info.shape[0], info.shape[1]);

            // create a vector for the data
            vector<double> vec(info.shape[0]*info.shape[1]);

            // copy the data over
            memcpy(vec.data(), reinterpret_cast<double *>(info.ptr),
                info.size * sizeof(double));

            // load in data
            m.load_raw(vec);

            // return the matrix
            return m;
        }),
            "Initializes a matrix from a buffer like a numpy array.\n\n"
            ":param b buffer: buffer to create a matrix from",
            "b"_a)
        .def("get_element", py::overload_cast<int, int>(&matrix_2::get_element),
            "Gets an element at a specific (i,j) location in the matrix.\n\n"
            ":param int i: location for the first index\n:param int j: location"
            " for the second index\n:return: returns the element at location "
            "(i, j) in the matrix\n:rtype: float",
            "i"_a, "j"_a)
        .def("get_element", py::overload_cast<int, int>(&matrix_2::get_element,
            py::const_),
            "Gets an element at a specific (i,j) location in the matrix.\n\n"
            ":param int i: location for the first index\n:param int j: location"
            " for the second index\n:return: returns the element at location "
            "(i, j) in the matrix\n:rtype: float (const)",
            "i"_a, "j"_a)
        .def("set_element", &matrix_2::set_element,
            "Sets an element at a specifi (i, j) location in the matrix.\n\n"
            ":param int i: location for the first index\n:param int j: location"
            " for the second index")
        .def("get_matrix", &matrix_2::get_matrix,
            "Get the internal matrix.\n\n:return: returns the internally stored"
            " matrix\n:rtype: List[List[float]]")
        .def_buffer([](matrix_2 &m) -> py::buffer_info {
            
            // get sizes and data
            auto [row, col] = m.get_size();
            auto ptr_val = m.get_matrix().data();

            // return buffer from this info
            return py::buffer_info(
                ptr_val, sizeof(double),
                py::format_descriptor<double>::format(), 2, {row, col},
                {sizeof(double) * col, sizeof(double)}
            );
        });

    // define the geometry class submodule
    auto geo = m.def_submodule("geometries", "Geometry classes for the lagmodel"
        " module");

    // define the base class type
    typedef Geometry<vector, double> geo_base;

    // define the geometry base class with it's methods
    py::class_<geo_base>(geo, "Geometry",
            "Defines a Geometry base class")
        .def(py::init<vector<double>, vector<double>, double>(),
            "Initializes the geometry from the geometry parameters.\n\n:param "
            "List[float] r: radial bins\n:param List[float] r_area: radial area"
            " bins\n:param float r_cor: radius of the corona",
            "r"_a, "r_area"_a, "r_cor"_a)
        .def("get_rcor", &geo_base::get_rcor,
            "Gets the radius of the corona.\n\n:return: returns the radius of "
            "the corona\n:rtype: float")
        .def("get_omega_cor", &geo_base::get_omega_cor,
            "Gets the solid angle array.\n\n:return: returns the solid angle\n"
            ":rtype: List[float]")
        .def("get_frad_disktocor", &geo_base::get_frad_disktocor,
            "Gets the disk to corona fraction array.\n\n:return: returns the "
            "disk to corona fraction\n:rtype: List[float]")
        .def("get_frad_cortodisk", &geo_base::get_frad_cortodisk,
            "Gets the corona to disk fraction array.\n\n:return: returns the "
            "corona to disk fraction\n:rtype: List[float]")
        .def("get_omega_cor_area", &geo_base::get_omega_cor_area,
            "Gets the solid angle multiplied with the radial area bins array."
            "\n\n:return: returns the solid angle multiplied with the radial "
            "area bins\n:rtype: List[float]")
        .def("get_frad_disktocor_area", &geo_base::get_frad_disktocor_area,
            "Gets the disk to corona fraction multiplied with the radial area "
            "bins array.\n\n:return: returns the disk to corona fraction "
            "multiplied with the radial area bins\n:rtype: List[float]")
        .def("get_frad_cortodisk_area", &geo_base::get_frad_cortodisk_area,
            "Gets the corona to disk fraction multiplied with the radial area "
            "bins array.\n\n:return: returns the corona to disk fraction "
            "multiplied with the radial area bins\n:rtype: List[float]");

    // define the an_sphere class type
    typedef AN_Sphere<vector, double> geo_ansphere;

    // define the an_sphere class
    py::class_<geo_ansphere, geo_base>(geo, "AN_Sphere",
            "Defines a Analytical Sphere class")
        .def(py::init<vector<double>, vector<double>, double>(),
            "Initializes the geometry from the geometry parameters.\n\n:param "
            "List[float] r: radial bins\n:param List[float] r_area: radial area"
            " bins\n:param float r_cor: radius of the corona",
            "r"_a, "r_area"_a, "r_cor"_a);

    // define the bknpow_emiss class type
    typedef BKNpow_Emiss<vector, double> geo_bknpowemiss;

    // define the bknpow_emiss class
    py::class_<geo_bknpowemiss, geo_base>(geo, "BKNpow_Emiss",
            "Defines a Broken Power-law Emission class")
        .def(py::init<vector<double>, vector<double>, double, double, double,
            double, double, double, double>(),
            "Initializes the geometry from the geometry parameters.\n\n:param "
            "List[float] r: radial bins\n:param List[float] r_area: radial area"
            " bins\n:param float r_cor: radius of the corona\n:param float "
            "r_bk1: radius of the break point 1\n:param float r_bk2: radius of "
            "the break point 2\n:param float em_ind1: induced emission of the "
            "first section\n:param float em_ind2: induced emission of the "
            "second section\n:param float em_ind3: induced emission of the "
            "third section\n:param float total_cortodisk: total corona to disk "
            "fraction",
            "r"_a, "r_area"_a, "r_cor"_a, "r_bk1"_a, "r_bk2"_a, "em_ind1"_a,
            "em_ind2"_a, "em_ind3"_a, "total_cortodisk"_a);

    // define the cylinder class type
    typedef Cylinder<vector, double> geo_cylinder;

    // define the cylinder class
    py::class_<geo_cylinder, geo_base>(geo, "Cylinder",
            "Defines a Cylinder class")
        .def(py::init<vector<double>, vector<double>, double, double, int,
            int>(),
            "Initializes the geometry from the geometry parameters.\n\n:param "
            "List[float] r: radial bins\n:param List[float] r_area: radial area"
            " bins\n:param float r_cor: radius of the corona\n:param float "
            "h_cor: height of the corona\n:param int nz: number of grid "
            "elements in the z direction\n:param int nphi: number of grid "
            "elements in the phi direction",
            "r"_a, "r_area"_a, "r_cor"_a, "h_cor"_a, "nz"_a, "nphi"_a)
        .def(py::init<vector<double>, vector<double>, double, double, double,
            double>(),
            "Initializes the geometry from the geometry parameters.\n\n:param "
            "List[float] r: radial bins\n:param List[float] r_area: radial area"
            " bins\n:param float r_cor: radius of the corona\n:param float "
            "h_cor: height of the corona\n:param float nz: number of grid "
            "elements in the z direction\n:param float nphi: number of grid "
            "elements in the phi direction",
            "r"_a, "r_area"_a, "r_cor"_a, "h_cor"_a, "nz"_a, "nphi"_a);

    // define the inv_cone class type
    typedef Inv_Cone<vector, double> geo_invcone;

    // define the inv_cone class
    py::class_<geo_invcone, geo_base>(geo, "Inv_Cone",
            "Defines a Inv_Cone class")
        .def(py::init<vector<double>, vector<double>, double, double, double,
            int, int>(),
            "Initializes the geometry from the geometry parameters.\n\n:param "
            "List[float] r: radial bins\n:param List[float] r_area: radial area"
            " bins\n:param float r_cor: radius of the corona\n:param float "
            "h_cor: height of the corona\n:param float r_top: radius at the top"
            "\n:param int nphi: number of grid elements in the phi direction\n"
            ":param int nz: number of grid elements in the z direction",
            "r"_a, "r_area"_a, "r_cor"_a, "h_cor"_a, "r_top"_a, "nphi"_a, "nz"_a)
        .def(py::init<vector<double>, vector<double>, double, double, double,
            double, double>());

    // define the sphere class type
    typedef Sphere<vector, double> geo_sphere;

    // define the sphere class
    py::class_<geo_sphere, geo_base>(geo, "Sphere",
            "Defines a Sphere class")
        .def(py::init<vector<double>, vector<double>, double, int, int>(),
            "Initializes the geometry from the geometry parameters.\n\n:param "
            "List[float] r: radial bins\n:param List[float] r_area: radial area"
            " bins\n:param float r_cor: radius of the corona\n:param int nphi: "
            "number of grid elements in the phi direction\n:param int ntheta: "
            "number of grid elements in the theta direction",
            "r"_a, "r_area"_a, "r_cor"_a, "nphi"_a, "nz"_a)
        .def(py::init<vector<double>, vector<double>, double, double, double>(),
            "Initializes the geometry from the geometry parameters.\n\n:param "
            "List[float] r: radial bins\n:param List[float] r_area: radial area"
            " bins\n:param float r_cor: radius of the corona\n:param float "
            "nphi: number of grid elements in the phi direction\n:param float "
            "ntheta: number of grid elements in the theta direction",
            "r"_a, "r_area"_a, "r_cor"_a, "nphi"_a, "nz"_a);
}