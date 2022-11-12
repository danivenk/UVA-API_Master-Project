#include "functions.h"
#include <fstream>
#include <json.hpp>
#include <ctime>
#include <cstdlib>
#include <sstream>

using json = nlohmann::json;

template <class T>
ostream& operator<<(ostream& os, const list<T> list);
template <class T>
ostream& operator<<(ostream& os, const Matrix<list, T> m);
int myrandom (int i);
template <class T>
json to_json(complex<T> c);
template <class T>
list<json> to_json(list<complex<T>> c_list);
template <class T>
string to_string(complex<T> c);
template <class T>
list<string> to_string(list<complex<T>> c_list);


void find_nearest(json &J, string f_name);
void calc_dispfrac(json &J, string f_name);
void calc_illumination_fracs(json &J, string f_name);
void calc_timing_params(json &J, string f_name);
void calc_propagation_params(json &J, string f_name);
void calc_radial_time_response(json &J, string f_name);
void calc_irfs_mono(json &J, string f_name);
void linear_rebin_irf(json &J, string f_name);
void calc_cross_psd(json &J, string f_name);
void lorentz_qold(json &J, string f_name);
void lorentz_q(json &J, string f_name);
void calculate_stprod_mono(json &J, string f_name);

int main() {
    // time(0)
    srand (unsigned (0));
    list<string> functions = {"find_nearest", "calc_dispfrac",
        "calc_illumination_fracs", "calc_timing_params",
        "calc_propagation_params", "calc_radial_time_response",
        "calc_irfs_mono", "linear_rebin_irf", "calc_cross_psd", "lorentz_qold",
        "lorentz_q", "calculate_stprod_mono"};
    typename list<string>::iterator function = functions.begin();

    json J;

    find_nearest(J, *(function++));
    calc_dispfrac(J, *(function++));
    calc_illumination_fracs(J, *(function++));
    calc_timing_params(J, *(function++));
    calc_propagation_params(J, *(function++));
    calc_radial_time_response(J, *(function++));
    calc_irfs_mono(J, *(function++));
    linear_rebin_irf(J, *(function++));
    calc_cross_psd(J, *(function++));
    lorentz_qold(J, *(function++));
    lorentz_q(J, *(function++));
    calculate_stprod_mono(J, *(function++));

    // cout << "next function: " << *(function) << endl;

    ofstream out("cpp_out.txt");
    out << J << endl;
    out.close();

    return 0;
}

void find_nearest(json &J, string f_name) {
    list<double> array;
    double value = .36;

    for (int i = 0; i < 10; i++) {
        array.push_back((double) i/10);
    }

    auto [found, idx] = find_nearest<double>(array, value);
    J[f_name]["test_1"] = {{"input", {array, value}}, {"output", {found, idx}}};

    value = -3.4;

    tie(found, idx) = find_nearest<double>(array, value);
    J[f_name]["test_2"] = {{"input", {array, value}}, {"output", {found, idx}}};

    value = 1.4;

    tie(found, idx) = find_nearest<double>(array, value);
    J[f_name]["test_3"] = {{"input", {array, value}}, {"output", {found, idx}}};

    value = 1e100;

    tie(found, idx) = find_nearest<double>(array, value);
    J[f_name]["test_4"] = {{"input", {array, value}}, {"output", {found, idx}}};

    value = -1e100;

    tie(found, idx) = find_nearest<double>(array, value);
    J[f_name]["test_5"] = {{"input", {array, value}}, {"output", {found, idx}}};

    value = .7;

    tie(found, idx) = find_nearest<double>(array, value);
    J[f_name]["test_6"] = {{"input", {array, value}}, {"output", {found, idx}}};
}

void calc_dispfrac(json &J, string f_name) {
    list<double> array_1, array_2;
    double value_1 = .37;
    double value_2 = 6.7;

    for (int i = 1; i <= 10; i++) {
        array_1.push_back((double) i/10);
        array_2.push_front((double) i/10);
    }

    auto [dispfrac, seed_frac_flow, heat_frac_flow] =
        calc_dispfrac<double>(array_1, array_2, value_1, value_1, value_2,
                              value_1, value_1, value_1);
    J[f_name]["test_1"] = {{"input", {array_1, array_2, value_1, value_1,
        value_2, value_1, value_1, value_1}}, {"output", {dispfrac,
        seed_frac_flow, heat_frac_flow}}};

    value_2 = 0;

    tie(dispfrac, seed_frac_flow, heat_frac_flow) =
        calc_dispfrac<double>(array_1, array_2, value_1, value_1, value_2,
                              value_1, value_1, value_1);
    J[f_name]["test_2"] = {{"input", {array_1, array_2, value_1, value_1,
        value_2, value_1, value_1, value_1}}, {"output", {dispfrac,
        seed_frac_flow, heat_frac_flow}}};

    value_2 = -4.9;

    tie(dispfrac, seed_frac_flow, heat_frac_flow) =
        calc_dispfrac<double>(array_1, array_2, value_1, value_1, value_2,
                              value_1, value_1, value_1);
    J[f_name]["test_3"] = {{"input", {array_1, array_2, value_1, value_1,
        value_2, value_1, value_1, value_1}}, {"output", {dispfrac,
        seed_frac_flow, heat_frac_flow}}};
}

void calc_illumination_fracs(json &J, string f_name) {
    list<double> array_1, array_2;
    double value_1 = .37;
    double value_2 = 6.7;
    double value_3 = .7;
    double value_4 = .13;
    typedef tuple<double> parm1;

    for (int i = 1; i <= 10; i++) {
        array_1.push_back((double) i/10);
        array_2.push_front((double) i/10);
    }

    auto r1 = make_tuple(array_1, array_2);

    parm1 parms1(value_1);
    auto f1 = tuple_cat(r1, parms1);

    auto geo1 = make_from_tuple<Geometry<list, double>>(f1);

    auto [omega_cor, frad_disktocor, frad_cortodisk] =
        calc_illumination_fracs<double>(geo1);
    J[f_name]["test_1"] = {{"input", {array_1, array_2, "geometry", parms1}},
        {"output", {omega_cor, frad_disktocor, frad_cortodisk}}};
    
    auto geo2 = make_from_tuple<AN_Sphere<list, double>>(f1);

    tie(omega_cor, frad_disktocor, frad_cortodisk) = 
        calc_illumination_fracs<double>(geo2);
    J[f_name]["test_2"] = {{"input", {array_1, array_2, "an_sphere", parms1}},
        {"output", {omega_cor, frad_disktocor, frad_cortodisk}}};

    typedef tuple<double, double, double, double, double, double, double> parm2;
    parm2 parms2(value_4, value_1, value_3, value_2/3, value_2/7, value_2/11,
        value_2);
    auto f2 = tuple_cat(r1, parms2);

    auto geo3 = make_from_tuple<BKNpow_Emiss<list, double>>(f2);

    tie(omega_cor, frad_disktocor, frad_cortodisk) = 
        calc_illumination_fracs<double>(geo3);
    J[f_name]["test_3"] = {{"input", {array_1, array_2, "bknpow_emiss", parms2}},
        {"output", {omega_cor, frad_disktocor, frad_cortodisk}}};

    parm2 parms3(value_4, value_3, value_1, value_2/3, value_2/7, value_2/11,
        value_2);
    auto f3 = tuple_cat(r1, parms3);

    auto geo4 = make_from_tuple<BKNpow_Emiss<list, double>>(f3);

    tie(omega_cor, frad_disktocor, frad_cortodisk) = 
        calc_illumination_fracs<double>(geo4);
    J[f_name]["test_4"] = {{"input", {array_1, array_2, "bknpow_emiss", parms3}},
        {"output", {omega_cor, frad_disktocor, frad_cortodisk}}};

    typedef tuple<double, int, int> parm3;
    parm3 parms4(value_1, (int) value_2, (int) value_2);
    auto f4 = tuple_cat(r1, parms4);

    auto geo5 = make_from_tuple<Sphere<list, double>>(f4);

    tie(omega_cor, frad_disktocor, frad_cortodisk) = 
        calc_illumination_fracs<double>(geo5);
    J[f_name]["test_5"] = {{"input", {array_1, array_2, "sphere", parms4}},
        {"output", {omega_cor, frad_disktocor, frad_cortodisk}}};

    typedef tuple<double, double, int, int> parm4;
    parm4 parms5(value_1, 3*value_1, (int) value_2, (int) value_2);
    auto f5 = tuple_cat(r1, parms5);

    auto geo6 = make_from_tuple<Cylinder<list, double>>(f5);

    tie(omega_cor, frad_disktocor, frad_cortodisk) = 
        calc_illumination_fracs<double>(geo6);
    J[f_name]["test_6"] = {{"input", {array_1, array_2, "cylinder", parms5}},
        {"output", {omega_cor, frad_disktocor, frad_cortodisk}}};

    typedef tuple<double, double, double, int, int> parm5;
    parm5 parms6(value_4, value_1, value_3, (int) value_2, (int) value_2);
    auto f6 = tuple_cat(r1, parms6);

    auto geo7 = make_from_tuple<Inv_Cone<list, double>>(f6);

    tie(omega_cor, frad_disktocor, frad_cortodisk) = 
        calc_illumination_fracs<double>(geo7);
    J[f_name]["test_7"] = {{"input", {array_1, array_2, "inv_cone", parms6}},
        {"output", {omega_cor, frad_disktocor, frad_cortodisk}}};
}

void calc_timing_params(json &J, string f_name) {
    list<double> array_1, array_2;
    double value_1 = .37;
    double value_2 = 6.7;
    double value_3 = .7;
    double value_4 = .13;
    int value_5 = 1;
    int value_6 = 4;

    for (int i = 1; i <= 10; i++) {
        array_1.push_back((double) i/10);
    }

    list<double> parms = {value_1};
    tuple<double, double> tparms1(value_2/1.3, value_3/1.4);
    tuple<double, double> tparms2(value_2/1.4, value_4/1.3);

    list<double>::iterator p = parms.begin();

    auto [tau, lfreq, q, rms] =
        calc_timing_params<double>(array_1, value_5, value_1, value_6, value_4,
            tparms1, tparms2, "", parms);
    J[f_name]["test_01"] = {{"input", {array_1, value_5, value_1, value_6,
        value_4, tparms1, tparms2, "", parms}}, {"output", {tau, lfreq, q,
        rms}}};

    parms.push_back(value_2);
    parms.push_back(value_2/1.7);

    tie(tau, lfreq, q, rms) =
        calc_timing_params<double>(array_1, 7 * value_5, value_1, value_6,
            value_4, tparms1, tparms2, "continuous", parms);
    J[f_name]["test_02"] = {{"input", {array_1, 7 * value_5, value_1, value_6,
        value_4, tparms1, tparms2, "continuous", parms}}, {"output", {tau,
        lfreq, q, rms}}};

    *next(p, 1) *= -1;

    tie(tau, lfreq, q, rms) =
        calc_timing_params<double>(array_1, 7 * value_5, value_1, value_6,
            value_4, tparms1, tparms2, "continuous", parms);
    J[f_name]["test_03"] = {{"input", {array_1, 7 * value_5, value_1, value_6,
        value_4, tparms1, tparms2, "continuous", parms}}, {"output", {tau,
        lfreq, q, rms}}};

    *next(p, 1) *= -1;
    *next(p, 2) *= -1;

    tie(tau, lfreq, q, rms) =
        calc_timing_params<double>(array_1, 7 * value_5, value_1, value_6,
            value_4, tparms1, tparms2, "continuous", parms);
    J[f_name]["test_04"] = {{"input", {array_1, 7 * value_5, value_1, value_6,
        value_4, tparms1, tparms2, "continuous", parms}}, {"output", {tau,
        lfreq, q, rms}}};

    *next(p, 1) *= -1;

    tie(tau, lfreq, q, rms) =
        calc_timing_params<double>(array_1, 7 * value_5, value_1, value_6,
            value_4, tparms1, tparms2, "continuous", parms);
    J[f_name]["test_05"] = {{"input", {array_1, 7 * value_5, value_1, value_6,
        value_4, tparms1, tparms2, "continuous", parms}}, {"output", {tau,
        lfreq, q, rms}}};

    *next(p, 1) *= -1;
    *next(p, 2) *= -1;

    tie(tau, lfreq, q, rms) =
        calc_timing_params<double>(array_1, 7 * value_5, value_1, value_6,
            value_4, tparms1, tparms2, "multi_frequency", parms);
    J[f_name]["test_06"] = {{"input", {array_1, 7 * value_5, value_1, value_6,
        value_4, tparms1, tparms2, "multi_frequency", parms}}, {"output", {tau,
        lfreq, q, rms}}};

    *p *= 3;

    tie(tau, lfreq, q, rms) =
        calc_timing_params<double>(array_1, 7 * value_5, value_1, value_6,
            value_4, tparms1, tparms2, "multi_frequency", parms);
    J[f_name]["test_07"] = {{"input", {array_1, 7 * value_5, value_1, value_6,
        value_4, tparms1, tparms2, "multi_frequency", parms}}, {"output", {tau,
        lfreq, q, rms}}};

    *p = 5 * value_1;

    tie(tau, lfreq, q, rms) =
        calc_timing_params<double>(array_1, 7 * value_5, value_1, value_6,
            value_4, tparms1, tparms2, "multi_frequency", parms);
    J[f_name]["test_08"] = {{"input", {array_1, 7 * value_5, value_1, value_6,
        value_4, tparms1, tparms2, "multi_frequency", parms}}, {"output", {tau,
        lfreq, q, rms}}};

    *p = value_1;

    list<double> parms2(parms.begin(), parms.end());

    *next(parms2.begin(), 0) *= 3;

    parms2.insert(parms2.end(), parms2.begin(), parms2.end());
    parms2.insert(parms2.end(), parms2.begin(), next(parms2.begin(), 3));

    *next(parms2.begin(), 6) = 5 * value_1;

    tie(tau, lfreq, q, rms) =
        calc_timing_params<double>(array_1, 7 * value_5, value_1, value_6,
            value_4, tparms1, tparms2, "multi_frequency", parms2);
    J[f_name]["test_09"] = {{"input", {array_1, 7 * value_5, value_1, value_6,
        value_4, tparms1, tparms2, "multi_frequency", parms2}}, {"output", {tau,
        lfreq, q, rms}}};

    tie(tau, lfreq, q, rms) =
        calc_timing_params<double>(array_1, 7 * value_5, value_1, value_6,
            value_4, tparms1, tparms2, "multi_radius", parms);
    J[f_name]["test_10"] = {{"input", {array_1, 7 * value_5, value_1, value_6,
        value_4, tparms1, tparms2, "multi_radius", parms}}, {"output", {tau,
        lfreq, q, rms}}};

    *p = value_3 * 1.3;

    tie(tau, lfreq, q, rms) =
        calc_timing_params<double>(array_1, 7 * value_5, value_1, value_6,
            value_4, tparms1, tparms2, "multi_radius", parms);
    J[f_name]["test_11"] = {{"input", {array_1, 7 * value_5, value_1, value_6,
        value_4, tparms1, tparms2, "multi_radius", parms}}, {"output", {tau,
        lfreq, q, rms}}};

    for (int i = 0; i < 2; i++) {
        *next(parms2.begin(), 3*i) = value_3 * 1.3;
    }
    *next(parms2.begin(), 6) = value_3 * 1.2;

    tie(tau, lfreq, q, rms) =
        calc_timing_params<double>(array_1, 7 * value_5, value_1, value_6,
            value_4, tparms1, tparms2, "multi_radius", parms2);
    J[f_name]["test_12"] = {{"input", {array_1, 7 * value_5, value_1, value_6,
        value_4, tparms1, tparms2, "multi_radius", parms2}}, {"output", {tau,
        lfreq, q, rms}}};

    parms.push_front(*p);
    p = parms.begin();
    *next(p, 1) = value_1;

    tie(tau, lfreq, q, rms) =
        calc_timing_params<double>(array_1, 7 * value_5, value_1, value_6,
            value_4, tparms1, tparms2, "multi_radius_frequency", parms);
    J[f_name]["test_13"] = {{"input", {array_1, 7 * value_5, value_1, value_6,
        value_4, tparms1, tparms2, "multi_radius_frequency", parms}},
        {"output", {tau, lfreq, q, rms}}};

    *p = value_3 * .9;

    tie(tau, lfreq, q, rms) =
        calc_timing_params<double>(array_1, 7 * value_5, value_1, value_6,
            value_4, tparms1, tparms2, "multi_radius_frequency", parms);
    J[f_name]["test_14"] = {{"input", {array_1, 7 * value_5, value_1, value_6,
        value_4, tparms1, tparms2, "multi_radius_frequency", parms}},
        {"output", {tau, lfreq, q, rms}}};

    for (int i = 0; i < 2; i++) {
        parms2.insert(next(parms2.begin(), 3*i + i), value_1);
    }
    
    parms2.insert(next(parms2.begin(), 8), value_3 * .9);

    tie(tau, lfreq, q, rms) =
        calc_timing_params<double>(array_1, 7 * value_5, value_1, value_6,
            value_4, tparms1, tparms2, "multi_radius_frequency", parms2);
    J[f_name]["test_14"] = {{"input", {array_1, 7 * value_5, value_1, value_6,
        value_4, tparms1, tparms2, "multi_radius_frequency", parms2}},
        {"output", {tau, lfreq, q, rms}}};
}

void calc_propagation_params(json &J, string f_name) {
    list<double> array_1, array_2;
    double value_1 = .37;
    double value_2 = 6.7;
    double value_3 = .7;
    double value_4 = .13;

    for (int i = 1; i <= 10; i++) {
        if (i != 10) {
            array_1.push_back((double) i/10);
        }
        array_2.push_back((double) (i-.5)/10);
    }

    tuple<double, double> parms1(value_2/1.3, value_3/1.4);
    tuple<double, double> parms2(value_2/1.4, value_4/1.3);

    auto deltau = calc_propagation_params<double>(array_1, array_2, value_1,
        parms1, parms2);
    J[f_name]["test_1"] = {{"input", {array_1, array_2, value_1, parms1,
        parms2}}, {"output", {deltau}}};
}

void calc_radial_time_response(json &J, string f_name) {
    list<double> array_1, array_2, array_3, array_4, array_5, array_6;
    double value_1 = .37;
    double value_2 = 6.7;
    double value_3 = .7;
    double value_4 = .13;

    for (int i = 1; i < 10; i++) {
        array_1.push_back((double) i/10);
        array_2.push_back((double) (myrandom(10))/10);
        array_3.push_back((double) (myrandom(10))/10);
        array_4.push_back((double) (myrandom(10))/10);
        array_5.push_back((double) (myrandom(10))/10);
        array_6.push_back((double) (myrandom(10))/10);
    }

    auto [ldisk_disp, lseed_disp, lheat, ldisk_rev, lseed_rev] = 
        calc_radial_time_response<double>(array_1, value_1, array_2, array_3,
            array_4, array_5, array_6);
    J[f_name]["test_1"] = {{"input", {array_1, value_1, array_2, array_3,
        array_4, array_5, array_6}}, {"output", {ldisk_disp, lseed_disp, lheat,
        ldisk_rev, lseed_rev}}};
}

void calc_irfs_mono(json &J, string f_name) {
    list<double> array_1, array_2, array_3, array_4, array_5, array_6;
    double value_1 = .37;
    double value_2 = 6.7;
    double value_3 = .7;
    double value_4 = .13;

    for (int i = 1; i < 10; i++) {
        array_1.push_back((double) i/10);
        array_2.push_back((double) (myrandom(10))/10);
        array_3.push_back((double) (myrandom(10))/10);
        array_4.push_back((double) (myrandom(10))/10);
        array_5.push_back((double) (myrandom(10))/10);
        array_6.push_back((double) (myrandom(10))/10);
    }

    tuple<double, double> parms(value_2/1.3, value_3/1.4);

    auto [gamma_mean1, gamma_irf1, flux_irf1, disk_irf1, seed_irf1] = 
        calc_irfs_mono<double, Array<list, double>>(parms, value_3,
            array_1, array_2, array_3, array_4, array_5, array_6);
    J[f_name]["test_1"] = {{"input", {parms, value_3, array_1, array_2, array_3,
        array_4, array_5, array_6}}, {"output", {gamma_mean1, gamma_irf1,
        flux_irf1.get_matrix(), disk_irf1, seed_irf1}}};

    auto [gamma_mean2, gamma_irf2, flux_irf2, disk_irf2, seed_irf2] =
        calc_irfs_mono<double, Nested_Array<list, double>>(parms, value_3,
            array_1, array_2, array_3, array_4, array_5, array_6);
    J[f_name]["test_2"] = {{"input", {parms, value_3, array_1, array_2, array_3,
        array_4, array_5, array_6}}, {"output", {gamma_mean2, gamma_irf2,
        flux_irf2.get_matrix(), disk_irf2, seed_irf2}}};
}

void linear_rebin_irf(json &J, string f_name) {
    list<double> array_1, array_2, array_3, array_4, array_5, array_6;
    double value_1 = .37;
    double value_2 = 6.7;
    double value_3 = .7;
    double value_4 = .13;
    int value_5 = 1;
    int value_6 = 4;

    for (int i = 1; i < 10; i++) {
        array_1.push_back((double) i/10);
        array_2.push_back((double) (myrandom(10))/10);
        array_3.push_back((double) (myrandom(10))/10);
        array_4.push_back((double) (myrandom(20))/10 - 1);
        array_5.push_back((double) (myrandom(10))/10);
        array_6.push_back((double) (myrandom(10))/10);
    }

    for (int i = 0; i < 1; i++) {
        array_1.pop_back();
    }

    auto [rebinned_irf, flux_outer] = 
        linear_rebin_irf<double>(value_4, 7 * value_5, array_6, array_4,
            array_3, array_2, 7 * value_5);
    J[f_name]["test_1"] = {{"input", {value_4, 7 * value_5, array_6, array_4,
       array_3, array_2, 7 * value_5}}, {"output", {rebinned_irf, flux_outer}}};

    tie(rebinned_irf, flux_outer) = 
        linear_rebin_irf<double>(value_4, 7 * value_5, array_6, array_4,
            array_1, array_2, 7 * value_5);
    J[f_name]["test_2"] = {{"input", {value_4, 7 * value_5, array_6, array_4,
       array_1, array_2, 7 * value_5}}, {"output", {rebinned_irf, flux_outer}}};
}

void calc_cross_psd(json &J, string f_name) {
    list<double> array_1, array_2, array_3, array_4, array_5, array_6, array_7,
        array_8, array_9;
    double value_1 = .37;
    double value_2 = 6.7;
    double value_3 = .7;
    double value_4 = .13;
    int value_5 = 1;
    int value_6 = 4;

    for (int i = 1; i < 10; i++) {
        array_1.push_back((double) i/10);
        array_2.push_back((double) (myrandom(10))/10);
        array_3.push_back((double) (myrandom(10))/10);
        array_4.push_back((double) (myrandom(10))/10);
        array_5.push_back((double) (myrandom(10))/10);
        array_6.push_back((double) (myrandom(10))/10);
        array_7.push_back((double) (myrandom(10))/10);
        array_8.push_back((double) (myrandom(10))/10);
        double x = (double) (myrandom(10))/10;
        array_9.push_back(x == 0 ? x + .1 : x);
    }

    for (int i = 1; i < 10; i++) {
        array_2.push_back((double) (myrandom(10)/10));
        array_3.push_back((double) (myrandom(10)/10));
    }

    auto [wt_cross_spec, wt_pow_spec_ci, wt_pow_spec_ref, mod_sig_psd]
        = calc_cross_psd<double>(array_7, value_4, array_2, array_3, 
            array_5, array_6, 5 * value_5, array_8, array_9, array_1);
    
    list<json> wt_cross_spec_s = to_json(wt_cross_spec);

    J[f_name]["test_1"] = {{"input", {array_7, value_4, array_2, array_3,
        array_5, array_6, 5 * value_5, array_8, array_9, array_1}}, {"output",
        {wt_cross_spec_s, wt_pow_spec_ci, wt_pow_spec_ref, mod_sig_psd}}};
}

void lorentz_qold(json &J, string f_name) {
    list<double> array_1;
    double value_1 = .37;
    double value_2 = 6.7;
    double value_3 = .7;

    for (int i = 1; i < 10; i++) {
        array_1.push_back((double) i/10);
    }

    auto lorentz = lorentz_q<double>(array_1, value_1, value_2, value_3);

    J[f_name]["test_1"] = {{"input", {array_1, value_1, value_2, value_3}},
        {"output", {lorentz}}};
}

void lorentz_q(json &J, string f_name) {
    list<double> array_1, array_2, array_3, array_4;

    for (int i = 1; i < 10; i++) {
        array_1.push_back((double) i/10);
    }

    for (int i = 1; i < 10; i++) {
        array_2.push_back((double) (myrandom(10)));
        array_4.push_back((double) (myrandom(10)));
        double x = (double) (myrandom(10))/10;
        array_3.push_back(x == 0 ? x + .1 : x);
    }

    auto lorentz_0 = lorentz_q<double, Array<list, double>>(array_1,
        array_2, array_3, array_4);

    J[f_name]["test_1"] = {{"input", {array_1, array_2, array_3, array_4}},
        {"output", {lorentz_0.get_matrix()}}};

    auto lorentz_1 = lorentz_q<double, Nested_Array<list, double>>(array_1,
        array_2, array_3, array_4);

    J[f_name]["test_2"] = {{"input", {array_1, array_2, array_3, array_4}},
        {"output", {lorentz_1.get_matrix()}}};
}

void calculate_stprod_mono(json &J, string f_name) {
    list<double> array_1, array_2, array_3, array_4, array_5, array_6, array_7,
        array_8, array_9;
    double value_1 = .37;
    double value_2 = 6.7;
    double value_3 = .7;
    double value_4 = .13;
    int value_5 = 1;
    int value_6 = 4;

    for (int i = 1; i < 10; i++) {
        array_1.push_back((double) i/10);
        array_2.push_back((double) (myrandom(10))/10);
        array_4.push_back((double) (myrandom(10))/10);
        array_5.push_back((double) (myrandom(10))/10);
        array_6.push_back((double) (myrandom(10))/10);
        array_7.push_back((double) (myrandom(10))/10);
        array_8.push_back((double) (myrandom(10))/10);
        double x = (double) (myrandom(10))/10;
        array_3.push_back(x == 0 ? x + .1 : x);
    }

    list<double> array_two = {(double) (myrandom(10))/10,
        (double) (myrandom(10))/10};

    Array<list, int> matrix_1(10, 2);
    Array<list, double> matrix_2(array_5, array_7);

    auto [x, y] = matrix_1.get_size();

    for (int i = 0; i < x; i++) {
        for (int j = 0; j < y; j++) {
            matrix_1.get_element(i, j) = myrandom(array_2.size());
        }
    }

    auto [freq1, phlag1, tlag1, psd_ci1, psd_ref1, mod_sig_psd1, irf_nbins1,
        irf_binedgefrac1, deltau_scale1, dt1, nirf1, ci_irf1, ci_mean1,
        ci_outer1] =
            calculate_stprod_mono<double, Array<list, double>,
                Array<list, int>>(value_6, array_2, matrix_1,
                    matrix_2, array_4, array_5, array_6, value_3, 7*value_5,
                    array_7, array_3, array_1, 2*value_4);

    J[f_name]["test_1"] = {{"input", {value_6, array_2, matrix_1.get_matrix(),
        matrix_2.get_matrix(), array_4, array_5, array_6, value_3, 7*value_5,
        array_7, array_3, array_1, 2*value_4}}, {"output", {freq1,
        phlag1.get_matrix(), tlag1.get_matrix(), psd_ci1.get_matrix(),
        psd_ref1.get_matrix(), mod_sig_psd1, irf_nbins1, irf_binedgefrac1,
        deltau_scale1, dt1, nirf1, ci_irf1.get_matrix(), ci_mean1, ci_outer1}}};

    Nested_Array<list, int> matrix_3(10, 2);
    Nested_Array<list, double> matrix_4(array_5, array_7);

    tie(x, y) = matrix_3.get_size();

    for (int i = 0; i < x; i++) {
        for (int j = 0; j < y; j++) {
            matrix_3.get_element(i, j) = myrandom(array_2.size());
        }
    }

    auto [freq2, phlag2, tlag2, psd_ci2, psd_ref2, mod_sig_psd2, irf_nbins2,
        irf_binedgefrac2, deltau_scale2, dt2, nirf2, ci_irf2, ci_mean2,
        ci_outer2] =
            calculate_stprod_mono<double, Nested_Array<list, double>,
                Nested_Array<list, int>>(value_6, array_2, matrix_3,
                    matrix_4, array_4, array_5, array_6, value_3, 7*value_5,
                    array_7, array_3, array_1, 2*value_4);

    J[f_name]["test_2"] = {{"input", {value_6, array_2, matrix_3.get_matrix(),
        matrix_4.get_matrix(), array_4, array_5, array_6, value_3, 7*value_5,
        array_7, array_3, array_1, 2*value_4}}, {"output", {freq2,
        phlag2.get_matrix(), tlag2.get_matrix(), psd_ci2.get_matrix(),
        psd_ref2.get_matrix(), mod_sig_psd2, irf_nbins2, irf_binedgefrac2,
        deltau_scale2, dt2, nirf2, ci_irf2.get_matrix(), ci_mean2, ci_outer2}}};
}

int myrandom (int i) { return std::rand()%i;}

template <class T>
json to_json(complex<T> c) {
    return {{"real", c.real()}, {"imag", c.imag()}};
}

template <class T>
list<json> to_json(list<complex<T>> c_list) {
    typename list<complex<T>>::iterator c;
    list<json> out;

    for (c = c_list.begin(); c != c_list.end(); c++) {
        out.push_back(to_json<T>(*c));
    }

    assert(c_list.size() == out.size());

    return out;
}

template <class T>
string to_string(complex<T> c) {
    stringstream ss;
    ss << c.real();
    if (c.imag() >= 0) {
        ss << "+";
    }
    ss << c.imag() << "j";
    return ss.str();
}

template <class T>
list<string> to_string(list<complex<T>> c_list) {
    typename list<complex<T>>::iterator c;
    list<string> out;
    
    for (c = c_list.begin(); c != c_list.end(); c++) {
        out.push_back(to_string<T>(*c));
    }
    
    assert(c_list.size() == out.size());
    
    return out;
}

template <class T>
ostream& operator<<(ostream& os, const list<T> list) {
    for (auto& i: list) {
        os << i << " ";
    }

    return os;
}

template <class T>
ostream& operator<<(ostream& os, const Matrix<list, T> m) {
    auto [x, y] = m.get_size();

    for (int i = 0; i < x; i++) {
        for (int j = 0; j < y; j++) {
            os << m.get_element(i, j) << " ";
        }
        os << endl;
    }

    return os;
}