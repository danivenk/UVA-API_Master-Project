#include "functions.h"
#include <fstream>
#include <json.hpp>

using json = nlohmann::json;

template <class T>
ostream& operator<<(ostream& os, const list<T> list);
void find_nearest(json &J, string f_name);
void calc_dispfrac(json &J, string f_name);
void calc_illumination_fracs(json &J, string f_name);

int main() {
    list<const char*> functions = {"find_nearest", "calc_dispfrac",
        "calc_illumination_fracs", "calc_timing_parms",
        "calc_propagation_parms", "calc_radial_time_reponse", "calc_irfs_mono",
        "linear_rebin_irf", "calc_cross_psd", "lorentz_qold", "lorentz_q",
        "calculate_stprod_mono"};
    typename list<const char*>::iterator function = functions.begin();

    json J;

    find_nearest(J, *(function++));
    calc_dispfrac(J, *(function++));

    cout << "next function: " << *(function) << endl;
    calc_illumination_fracs(J, *(function++));

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

template <class T>
ostream& operator<<(ostream& os, const list<T> list) {
    for (auto& i: list) {
        os << i << " ";
    }

    return os;
}