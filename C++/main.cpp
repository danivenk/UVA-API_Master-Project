#include <iostream>
#include <list>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <cassert>
#include <iterator>
#include <typeinfo>

#include "includes.h"
#include "geometries/bknpow_emiss.h"
#include "geometries/an_sphere.h"
#include "geometries/sphere.h"
#include "geometries/cylinder.h"
#include "geometries/inv_cone.h"

template <class T>
tuple<T, int> find_nearest(list<T> array, T value);
template <class T>
ostream& operator<<(ostream& os, const list<T> list);

int main() {

    srand(time(nullptr));

    list<double> _list = {}, _area = {};

    for (int i = 0; i < 10; i++) {
        // _list.push_back(rand());
        _list.push_back(i);
        _area.push_front(i);
    }

    using test_bnk_pow = tuple<double, double, double, double, double, double, double>;
    test_bnk_pow test1 = make_tuple(1,2,3,4,5,6,7);

    cout << "control" << endl;

    Geometry<list<double>, test_bnk_pow> a(_list, _area, test1);

    cout << a.get_omega_cor() << endl;
    cout << a.get_frad_disktocor() << endl;
    cout << a.get_frad_cortodisk() << endl;
    
    cout << "BKN Power Emssion" << endl;

    BKNpow_Emiss<list<double>, test_bnk_pow> b(_list, _area, test1);

    cout << b.get_omega_cor() << endl;
    cout << b.get_frad_disktocor() << endl;
    cout << b.get_frad_cortodisk() << endl;

    using test_an_sphere = tuple<double>;
    test_an_sphere test2 = make_tuple(1);

    cout << "AN Sphere" << endl;

    AN_Sphere<list<double>, test_an_sphere> c(_list, _area, test2);

    cout << c.get_omega_cor() << endl;
    cout << c.get_frad_disktocor() << endl;
    cout << c.get_frad_cortodisk() << endl;

    using test_sphere = tuple<double, double, double>;
    test_sphere test3 = make_tuple(1,2,3);

    cout << "Sphere" << endl;

    Sphere<list<double>, test_sphere> d(_list, _area, test3);

    cout << d.get_omega_cor() << endl;
    cout << d.get_frad_disktocor() << endl;
    cout << d.get_frad_cortodisk() << endl;

    using test_cylinder = tuple<double, double, double, double>;
    test_cylinder test4 = make_tuple(1,2,3,4);

    cout << "Cylinder" << endl;

    Cylinder<list<double>, test_cylinder> e(_list, _area, test4);

    cout << e.get_omega_cor() << endl;
    cout << e.get_frad_disktocor() << endl;
    cout << e.get_frad_cortodisk() << endl;

    using test_inv_cone = tuple<double, double, double, double, double>;
    test_inv_cone test5 = make_tuple(1,2,3,4,5);

    cout << "Inverse Cone" << endl;

    Inv_Cone<list<double>, test_inv_cone> f(_list, _area, test5);

    cout << f.get_omega_cor() << endl;
    cout << f.get_frad_disktocor() << endl;
    cout << f.get_frad_cortodisk() << endl;

    // A a(2);

    // cout << "A" << endl;
    // cout << a.get_a() << endl;
    // cout << square(a) << endl;

    // B b(2);

    // cout << "B" << endl;
    // cout << b.get_a() << endl;
    // cout << square(b) << endl;

    return 0;
}

template <class T>
tuple<T, int> find_nearest(list<T> array, T value) {
    
    cout << "list: " << array << endl;
    cout << "value: " << value << endl;

    T smallest = numeric_limits<T>::max();
    T small_diff = numeric_limits<T>::max();

    int i = 0;

    for (T &item: array) {
        T diff = abs(item - value);

        cout << smallest << " - " << item << endl;
        cout << small_diff << " _ " << diff << endl;

        if (diff < small_diff) {
            smallest = item;
            small_diff = diff;
        }

        i++;
    }
    
    return make_tuple(smallest, i);
}

template <class T>
tuple<list<T>, list<T>, list<T>> calc_dispfrac(list<T> rad, list<T> rad_area,
        T rin, T rcor, T seedff_norm, T seedff_ind, T heatff_norm,
        T heatff_ind) {

    assert (rad.size() == rad_area.size());

    list<T> dispfrac (rad.size(), 0);
    list<T> seed_frac_flow = (rad.size(), 0);
    list<T> heat_frac_flow = (rad.size(), 0);

    typename list<T>::iterator disp, seed, heat, r, area;

    for (disp = dispfrac.begin(), seed = seed_frac_flow.begin(),
            heat = heat_frac_flow.begin(), r = rad.begin(),
            area = rad_area.begin(); disp != dispfrac.end(),
            seed != seed_frac_flow.end(), heat != heat_frac_flow.end(),
            r != rad.end(), area != rad_area.end(); disp++, seed++, heat++,
            r++, area++) {
        *disp = *area * pow(*r, -3.0) * (1 - sqrt(rin/(*r)));

        if (*r <= rcor) {
            *heat = pow(heatff_norm*(*r/rin), heatff_ind);

            if (seedff_norm >= 0) {
                *seed = pow(seedff_norm*(*r/rin), heatff_ind);
            } else {
                *seed = 1.0 - *heat;
            }
        }
    }

    return make_tuple(dispfrac, seed_frac_flow, heat_frac_flow);
}

template <class T>
tuple<list<T>, list<T>, list<T>> calc_illumination_fracs(list<T> rad,
        list<T> rad_area, Geometry<list<T>, T> geomod, list<T> parms) {

    assert(rad.size() == rad_area.size());

    list<T> omega_cor(geomod.get_omega_cor()),
        frad_disktocor(geomod.get_frad_disktocor()),
        frad_cortodisk (geomod.get_frad_cortodisk());

    typename list<T>::iterator fcd, area;

    for (fcd = geomod.get_frad_cortodisk().begin(), area = rad_area.begin();
            fcd != geomod.get_frad_cortodisk().end(), area != rad_area.end();
            fcd++, area++) {
        *fcd *= *area;
    }

    return make_tuple(omega_cor, frad_disktocor, frad_cortodisk);
}

template <class T, class U>
tuple<list<T>, list<T>, list<T>, list<T>> calc_timing_parms(list<T> rad,
        int i_rsigmax, T rcor, int i_rcor, T t_scale, tuple<T, T> disk_tau_par,
        tuple<T, T> cor_tau_par, const char* lor_model, tuple<U> lor_par) {
    
    list<T> tau, lfreq, q_list, rms;

    typename list<T>::iterator r;
    int i;

    for (r = rad.begin(), i = 0; r != rad.end(), i < rad.size(); r++, i++) {
        if (*r > rcor) {
            tau.push_back(pow(get<0>(disk_tau_par) * (*r/rcor),get<1>(disk_tau_par))*pow(*r, 1.5));
        } else {
            tau.push_back(pow(get<0>(cor_tau_par) * (*r/rcor),get<1>(cor_tau_par))*pow(*r, 1.5));
        }
        
        if (lor_model == "continous") {
            assert(tuple_size<tuple<U>>{} == 3);
            lfreq.push_back(1/tau.back()*t_scale);
            q_list.push_back(get<0>(lor_par));

            if (i > i_rsigmax) {
                rms.push_back(0);
            } else {
                if (get<1>(lor_par) >= 0 && get<2>(lor_par) >= 0) {
                    if (*r > rcor) {
                        rms.push_back(sqrt(pow(get<1>(lor_par), 2)/
                            (i_rsigmax+1-i_rcor)));
                    } else {
                        rms.push_back(sqrt(pow(get<2>(lor_par), 2)/(i_rcor)));
                    }
                } if else (get<1>(lor_par) >= 0 && get<2>(lor_par) < 0) {
                    rms.push_back(sqrt(pow(get<1>(lor_par), 2)/(i_rsigmax+1)));
                } if else (get<1>(lor_par) < 0 && get<2>(lor_par) >= 0) {
                    if (*r > rcor) {
                        rms.push_back(-get<1>(lor_par));
                    } else {
                        rms.push_back(sqrt(pow(get<2>(lor_par), 2)/i_rcor));
                    }
                } else {
                    if (*r > rcor) {
                        rms.push_back(-get<1>(lor_par));
                    } else {
                        rms.push_back(-get<2>(lor_par));
                    }
                }
            }
        } else {
            lfreq.push_back(0);
            q_list.push_back(0);
            rms.push_back(0);
        }
    }

    r = rad.begin();

    typename list<T>::iterator t = tau.begin(), freq = lfreq.begin(),
        q = q_list.begin(), rms_value = rms.begin(); 

    for (int j = 0; j < tuple_size<U>{}; j += 3) {
        if (lor_model == "multi_frequency") {
            double ltau = 1/(get<j>(lor_par)*t_scale);
            if (ltau <= *(t + tau.size() - 1) && ltau >= *(tau.begin())) {
                auto [ltau2, i] = find_nearest<T>(tau, ltau);
                *(freq + i) = 1/(*(t + i) * t_scale);
                *(q + i) = get<j+1>(lor_par);
                if (i < i_rsigmax) {
                    *(rms_value + i) = get<j+2>(lor_par);
                }
            }
        }
        else if (lor_model == "multi_radius") {
            if (get<j>(lor_par) <= *(r + rad.size() - 1) &&
                    get<j>(lor_par) >= *r) {
                auto [lrad, i] = find_nearest<T>(rad, get<j>(lor_par));
                *(freq + i) = 1/(*(t + 1) * t_scale);
                *(q + i) = get<j+1>(lor_par);
                if (i < i_rsigmax) {
                    *(rms_value + i) = get<j+2>(lor_par);
                }
            }
        } else if (lor_model == "multi_radius_frequency") {
            if (get<j>(lor_par) <= *(r + rad.size() - 1) &&
                    get<j>(lor_par) >= *r) {
                auto [lrad, i] = find_nearest<T>(rad, get<j>(lor_par));
                if (*(r + i) > 0) {
                    for (int i_next = 1; i_next < i; i_next++) {
                        if (*(rms_value + i - i_next) == 0 &&
                                *(rms_value + i) > 0) {
                            i = i - i_next;
                        }
                    }
                    *(freq + i) = get<j+1>(lor_par);
                    *(q + i) = get<j+2>(lor_par);
                    if (i < i_rsigmax) {
                        *(rms + i) = get<j+3>(lor_par);
                    }

                    cout << "Signal radius for " << get<j+1>(lor_par) <<
                        " Hz is " << *(r + i) << endl;
                    cout << "with rms" << *(rms_value + i)
                }
            }
        }
    }

    return make_tuple(tau, lfreq, q_list, rms);
}

template <class T>
ostream& operator<<(ostream& os, const list<T> list) {
    for (auto& i: list) {
        os << i << " ";
    }

    return os;
}