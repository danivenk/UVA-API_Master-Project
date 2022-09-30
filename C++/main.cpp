#include "functions.h"
#include <fstream>
#include <json.hpp>

using json = nlohmann::json;

template <class T>
ostream& operator<<(ostream& os, const list<T> list);
void find_nearest(json &J, string f_name);
void calc_dispfrac(json &J, string f_name);
void calc_illumination_fracs(json &J, string f_name);
void calc_timing_params(json &J, string f_name);

int main() {
    list<const char*> functions = {"find_nearest", "calc_dispfrac",
        "calc_illumination_fracs", "calc_timing_params",
        "calc_propagation_parms", "calc_radial_time_reponse", "calc_irfs_mono",
        "linear_rebin_irf", "calc_cross_psd", "lorentz_qold", "lorentz_q",
        "calculate_stprod_mono"};
    typename list<const char*>::iterator function = functions.begin();

    json J;

    find_nearest(J, *(function++));
    calc_dispfrac(J, *(function++));
    calc_illumination_fracs(J, *(function++));

    cout << "next function: " << *(function) << endl;
    calc_timing_params(J, *(function++));

    // cout << J << endl;

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

    parm1 parms1(value_1);

    Geometry<list<double>, parm1> geo1(array_1, array_2, parms1);

    auto [omega_cor, frad_disktocor, frad_cortodisk] =
        calc_illumination_fracs<double, parm1>(geo1);
    J[f_name]["test_1"] = {{"input", {array_1, array_2, "geometry", parms1}},
        {"output", {omega_cor, frad_disktocor, frad_cortodisk}}};
    
    AN_Sphere<list<double>, parm1> geo2(array_1, array_2, parms1);

    tie(omega_cor, frad_disktocor, frad_cortodisk) = 
        calc_illumination_fracs<double, parm1>(geo2);
    J[f_name]["test_2"] = {{"input", {array_1, array_2, "an_sphere", parms1}},
        {"output", {omega_cor, frad_disktocor, frad_cortodisk}}};

    typedef tuple<double, double, double, double, double, double, double> parm2;
    parm2 parms2(value_4, value_1, value_3, value_2/3, value_2/7, value_2/11,
        value_2);

    BKNpow_Emiss<list<double>, parm2> geo3(array_1, array_2, parms2);

    tie(omega_cor, frad_disktocor, frad_cortodisk) = 
        calc_illumination_fracs<double, parm2>(geo3);
    J[f_name]["test_3"] = {{"input", {array_1, array_2, "bknpow_emiss", parms2}},
        {"output", {omega_cor, frad_disktocor, frad_cortodisk}}};

    parm2 parms3(value_4, value_3, value_1, value_2/3, value_2/7, value_2/11,
        value_2);

    BKNpow_Emiss<list<double>, parm2> geo4(array_1, array_2, parms3);

    tie(omega_cor, frad_disktocor, frad_cortodisk) = 
        calc_illumination_fracs<double, parm2>(geo4);
    J[f_name]["test_4"] = {{"input", {array_1, array_2, "bknpow_emiss", parms3}},
        {"output", {omega_cor, frad_disktocor, frad_cortodisk}}};

    typedef tuple<double, int, int> parm3;
    parm3 parms4(value_1, (int) value_2, (int) value_2);

    Sphere<list<double>, parm3> geo5(array_1, array_2, parms4);

    tie(omega_cor, frad_disktocor, frad_cortodisk) = 
        calc_illumination_fracs<double, parm3>(geo5);
    J[f_name]["test_5"] = {{"input", {array_1, array_2, "sphere", parms4}},
        {"output", {omega_cor, frad_disktocor, frad_cortodisk}}};

    typedef tuple<double, double, int, int> parm4;
    parm4 parms5(value_1, 3*value_1, (int) value_2, (int) value_2);

    Cylinder<list<double>, parm4> geo6(array_1, array_2, parms5);

    tie(omega_cor, frad_disktocor, frad_cortodisk) = 
        calc_illumination_fracs<double, parm4>(geo6);
    J[f_name]["test_6"] = {{"input", {array_1, array_2, "cylinder", parms5}},
        {"output", {omega_cor, frad_disktocor, frad_cortodisk}}};

    typedef tuple<double, double, double, int, int> parm5;
    parm5 parms6(value_4, value_1, value_3, (int) value_2, (int) value_2);

    Inv_Cone<list<double>, parm5> geo7(array_1, array_2, parms6);

    tie(omega_cor, frad_disktocor, frad_cortodisk) = 
        calc_illumination_fracs<double, parm5>(geo7);
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
        // array_2.push_front((double) i/10);
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

template <class T>
ostream& operator<<(ostream& os, const list<T> list) {
    for (auto& i: list) {
        os << i << " ";
    }

    return os;
}