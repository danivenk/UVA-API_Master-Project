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
#include "array.h"

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

    list<double> a_list, b_list, c_list;

    for (int i = 0; i < 10; i++) {
        a_list.push_back(i);
        if (i % 4 == 0) {
            b_list.push_back(i);
        }
        if (i % 4 == 1) {
            c_list.push_back(i);
        }
    }

    Nested_Array<list<double>, list<list<double>>, double> test(a_list, b_list);
    Array<list<double>, double> tset(a_list, b_list);
    // Array<list<double>, double> teet(a_list, c_list);

    for (int i = 0; i < 5; i++) {
        cout << "-";
    }
    cout << endl;

    Nested_Array<list<double>, list<list<double>>, double> teet = test + test;
    Array<list<double>, double> tsst = tset + tset;

    // cout << test << endl;
    // cout << tset << endl;
    cout << teet << endl;
    cout << tsst << endl;

    // for (a_val = a_list.begin(), b_val = b_list.begin(); a_val != a_list.end(),
    //         b_val != b_list.end(); a_val++, b_val++) {
    //     cout << *a_val << " - " << (*b_val-*(b_val + 1)) << endl;
    // }

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
            tau.push_back(pow(get<0>(disk_tau_par) * (*r/rcor), 
                get<1>(disk_tau_par))*pow(*r, 1.5));
        } else {
            tau.push_back(pow(get<0>(cor_tau_par) * (*r/rcor),
                get<1>(cor_tau_par))*pow(*r, 1.5));
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
                } else if (get<1>(lor_par) >= 0 && get<2>(lor_par) < 0) {
                    rms.push_back(sqrt(pow(get<1>(lor_par), 2)/(i_rsigmax+1)));
                } else if (get<1>(lor_par) < 0 && get<2>(lor_par) >= 0) {
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
            if (ltau <= *next(t, tau.size() - 1) && ltau >= *(tau.begin())) {
                auto [ltau2, i] = find_nearest<T>(tau, ltau);
                *next(freq, i) = 1/(*next(t, i) * t_scale);
                *next(q, i) = get<(const int) j+1>(lor_par);
                if (i < i_rsigmax) {
                    *next(rms_value, i) = get<(const int) j+2>(lor_par);
                }
            }
        }
        else if (lor_model == "multi_radius") {
            if (get<j>(lor_par) <= *(r + rad.size() - 1) &&
                    get<j>(lor_par) >= *r) {
                auto [lrad, i] = find_nearest<T>(rad, get<j>(lor_par));
                *next(freq, i) = 1/(*next(t, 1) * t_scale);
                *next(q, i) = get<(const int) j+1>(lor_par);
                if (i < i_rsigmax) {
                    *next(rms_value, i) = get<(const int) j+2>(lor_par);
                }
            }
        } else if (lor_model == "multi_radius_frequency") {
            if (get<j>(lor_par) <= *(r + rad.size() - 1) &&
                    get<j>(lor_par) >= *r) {
                auto [lrad, i] = find_nearest<T>(rad, get<j>(lor_par));
                if (*next(r, i) > 0) {
                    for (int i_next = 1; i_next < i; i_next++) {
                        if (*next(rms_value, i - i_next) == 0 &&
                                *next(rms_value, i) > 0) {
                            i = i - i_next;
                        }
                    }
                    *next(freq, i) = get<(const int) j+1>(lor_par);
                    *next(q, i) = get<(const int) j+2>(lor_par);
                    if (i < i_rsigmax) {
                        *next(rms_value, i) = get<(const int) j+3>(lor_par);
                    }

                    cout << "Signal radius for " << 
                        get<(const int) j+1>(lor_par) << " Hz is " <<
                        *next(r, i) << endl << "with rms" <<
                        *next(rms_value, i);
                }
            }
        }
    }

    return make_tuple(tau, lfreq, q_list, rms);
}

template <class T>
list <T> cacl_propagation_parms(list<T> rad, list<T> rad_edge, T rcor,
        tuple<T, T> disk_prop_par, tuple<T, T> cor_prop_par) {
    list<T> deltau;

    typename list<T>::iterator r, edge;

    for (r = rad.begin(), edge = rad_edge.begin(); r = rad.end(),
            edge = rad_edge.end(); r++, edge++) {
        if (*r <= rcor) {
            deltau.append(get<0>(cor_prop_par) * (*next(edge, 1) - *edge) * 
                pow(*r, 1/2) * pow(*r/rcor, get<1>(cor_prop_par)));
        } else {
            deltau.append(get<0>(disk_prop_par) * (*next(edge, 1) - *edge) * 
                pow(*r, 1/2) * pow(*r/rcor, get<1>(disk_prop_par)));
        }
    }

    return deltau;
}

template <class T>
tuple<list<T>, list<T>, list<T>, list<T>, list<T>> calc_radial_time_reponse(
        list<T> rad, T rcor, list<T> disp_frac, list<T> disktocor_frac,
        list<T> cortodisk_frac, list<T> seed_frac_flow, list<T> heat_frac_flow) {
    
    list<T> ldisk_disp, lseed_disp, lheat, ldisk_rev, lseed_rev;

    typename list<T>::iterator r, disp, fdc, fcd, seed, heat;

    double f_rev = 0, f_return = 0;

    for (fcd = cortodisk_frac.begin(), fdc = disktocor_frac.begin();
            fcd != cortodisk_frac.end(), fdc = disktocor_frac.end(); fcd++,
            fdc++) {
        f_rev += *fcd; f_return += (*fcd) * (*fdc);
    }

    for (r = rad.begin(), disp = disp_frac.begin(),
            fdc = disktocor_frac.begin(), seed = seed_frac_flow.begin(),
            heat = heat_frac_flow.begin(); r != rad.end(),
            disp != disp_frac.end(), fdc = disktocor_frac.end(),
            seed = seed_frac_flow.end(), heat = heat_frac_flow.end(); r++,
            disp++, fdc++, seed++, heat++) {
        if (*r > rcor) {
            ldisk_disp.push_back((*disp) * (1 - *fdc));
        } else {
            ldisk_disp.push_back(0);
        }

        lseed_disp.push_back((*disp) * (*fdc + *seed) * (1 - f_rev));
        lheat.push_back((*disp) * (*heat) * (1 - f_rev));

        lseed_rev.push_back(f_return * (1 - f_rev) * (*disp) *
            (*fdc + *seed + *heat)/(1 - f_return));
        ldisk_rev.push_back((f_rev - f_return) * (*disp) *
            (*fdc + *seed + *heat)/(1 - f_return));
    }

    return make_tuple(ldisk_disp, lseed_disp, lheat, ldisk_rev, lseed_rev);
}

template <class T, class U>
tuple<T, list<T>, U, list<T>, list<T>> calc_irfs_mono(tuple<T, T> gamma_par,
        T e_seed, list<T> energy, list<T> ldisk_disp, list<T> lseed_disp,
        list<T> lheat, list<T> ldisk_rev, list <T> lseed_rev) {

    double lheat_sum, lseed_sum, gamma_mean, u, v, n;
    list<T> gamma_irf, disk_irf, seed_irf, seed_div_sum_irf, heat_div_sum_irf,
        seed_column, heat_column;

    typename list<T>::iterator heat, seed_disp, seed_rev, disk_disp, disk_rev,
        energy_value;

    lheat_sum = 0, lseed_sum = 0;

    for (heat = lheat.begin(), seed_disp = lseed_disp.begin(),
            seed_rev = lseed_rev.begin(); heat = lheat.end(),
            seed_disp = lseed_disp.end(), seed_rev = lseed_rev.end(); heat++,
            seed_disp++, seed_rev++) {
        lheat_sum += *heat; lseed_sum += *seed_disp + *seed_rev;
    }

    gamma_mean = get<0>(gamma_par) * pow(lseed_sum/lheat_sum, get<1>(gamma_par));
    u = gamma_mean * get<1>(gamma_par);
    v = gamma_mean - 1;

    for (energy_value = energy.begin(); energy_value = energy.end();
            energy_value++) {
        seed_column.push_back(1-u*log((*energy_value)/e_seed) + (u/v));
        heat_column.push_back(u*log((*energy_value)/e_seed) - (u/v));
    }

    for (seed_disp = lseed_disp.begin(), seed_rev = lseed_rev.begin(),
            heat = lheat.begin(), disk_disp = ldisk_disp.begin(),
            disk_rev = ldisk_rev.begin(); seed_disp = lseed_disp.end(),
            seed_rev = lseed_rev.end(), heat = lheat.end(),
            disk_disp = ldisk_disp.end(), disk_rev = ldisk_rev.end();
            seed_disp++, seed_rev++, heat++, disk_disp++, disk_rev++) {
        gamma_irf.push_back(gamma_mean * get<1>(gamma_par) *
            ((*seed_disp + *seed_rev)/lseed_sum - *heat/lheat_sum));
        disk_irf.push_back(*disk_disp + *disk_rev);
        seed_irf.push_back(*seed_disp + *seed_rev);
        seed_div_sum_irf.push_back((*seed_disp + *seed_rev)/lseed_sum);
        heat_div_sum_irf.push_back(*heat/lheat);
    }

    U firf_seed(seed_column, seed_div_sum_irf),
        firf_heat(heat_column, heat_div_sum_irf);

    U flux_irf = firf_seed + firf_heat;

    make_tuple(gamma_mean, gamma_irf, disk_irf, seed_irf);
}

template <class T>
tuple<list<T>, T> linear_rebin_irf(T dt, int i_rsigmax, list<T> irf_nbins,
        list<T> irf_binedgefrac, list<T> input_irf, list<T> deltau_scale,
        int nirf) {
    list<T> rebinned_irf(nirf, 0), input_irf2(input_irf.size(), 0);

    typename list<T>::iterator ds = deltau_scale.begin(),
        irf = input_irf.begin(), irf2 = input_irf2.begin(),
        nbins = irf_nbins.begin(), rebinned, bef = irf_binedgefrac.begin();
    
    for (int i = i_rsigmax; i > 0; i--) {
        if (*next(ds, i) > 0) {
            *next(irf2, i) = *next(irf, i) / *next(ds, i);
        }
    }
    
    T irfbin = 0;

    for (int j = i; j <= i_rsigmax; j++) {
        irfbin += *next(nbins, j);
    }

    int irfbin_start = static_cast<int>(irfbin) -(*next(nbins, i));
    int irfbin_stop = irfbin - 1;

    for (rebinned = rebinned_irf.begin(); rebinned = rebinned_irf.end();
            rebinned++) {
        *rebinned = *next(irf2, i);
    }

    rebinned = rebinned_irf.begin();

    for (int i = i_rsigmax; i > 0; i--) {
        if (*next(ds, i) > 0) {
            if (i < i_rsigmax && *next(bef, i) > 0) {
                *next(rebinned, irfbin_start) = *next(irf2, i) +
                    (*next(bef, i) * (*next(irf2, i+1) - *next(irf2, i)))
            }
            if (i > 0 && *next(bef, i - 1) <= 0) {
                *next(rebinned, irfbin_stop) = *next(irf2, i) +
                    (*next(bef, i-1) * (*next(irf2, i) - *next(irf2, i-1)))
            }
        } else {
            *next(rebinned, irfbin_stop+1) += *next(irf2, i);
        }
    }

    *next(rebinned, irfbin_stop+1) /= dt;

    T flux_outer = 0

    if (i_rsigmax < input_irf.size() - 1) {
        for (int i = i_rsigmax; i < input_irf.size(); i++) {
            flux_outer += *next(irf, i);
        }
    }

    make_tuple(rebinned_irf, flux_outer);
} 

template <class T>
ostream& operator<<(ostream& os, const list<T> list) {
    for (auto& i: list) {
        os << i << " ";
    }

    return os;
}