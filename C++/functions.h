#ifndef functions_H
#define functions_H

#include <vector>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <cassert>
#include <iterator>

#include "includes.h"
#include "geometries/bknpow_emiss.h"
#include "geometries/an_sphere.h"
#include "geometries/sphere.h"
#include "geometries/cylinder.h"
#include "geometries/inv_cone.h"
#include "FFT/FFT.h"
#include "array.h"

template <class T>
tuple<T, int> find_nearest(list<T> array, T value);
template <class T>
tuple<list<T>, list<T>, list<T>> calc_dispfrac(list<T> rad, list<T> rad_area,
    T rin, T rcor, T seedff_norm, T seedff_ind, T heatff_norm, T heatff_ind);
template <class T, class U>
tuple<list<T>, list<T>, list<T>> calc_illumination_fracs(Geometry<list<T>, U> geomod);
template <class T>
tuple<list<T>, list<T>, list<T>, list<T>> calc_timing_params(list<T> rad,
    int i_rsigmax, T rcor, int i_rcor, T t_scale, tuple<T, T> disk_tau_par,
    tuple<T, T> cor_tau_par, string lor_model, list<T> lor_par);
template <class T>
list<T> calc_propagation_parms(list<T> rad, list<T> rad_edge, T rcor,
    tuple<T, T> disk_prop_par, tuple<T, T> cor_prop_par);
template <class T>
tuple<list<T>, list<T>, list<T>, list<T>, list<T>> calc_radial_time_reponse(
    list<T> rad, T rcor, list<T> disp_frac, list<T> disktocor_frac,
    list<T> cortodisk_frac, list<T> seed_frac_flow, list<T> heat_frac_flow);
template <class T, class U>
tuple<T, list<T>, U, list<T>, list<T>> calc_irfs_mono(tuple<T, T> gamma_par,
    T e_seed, list<T> energy, list<T> ldisk_disp, list<T> lseed_disp,
    list<T> lheat, list<T> ldisk_rev, list <T> lseed_rev);
template <class T>
tuple<list<T>, T> linear_rebin_irf(T dt, int i_rsigmax, list<T> irf_nbins,
    list<T> irf_binedgefrac, list<T> input_irf, list<T> deltau_scale, int nirf);
template <class T>
tuple<list<complex<T>>, list<T>, list<T>, list<T>> calc_cross_psd(list<T> freq,
    T dt, list<T> ci_irf, list<T> ref_irf, list <T> irf_nbins,
    list<T> deltau_scale, int i_rsigmax, list<T> f_pk, list<T> q, list<T> rms);
template <class T, class U>
U lorentz_q(list<T> f, list<T> f_pk, list<T> q, list<T> rms);
template<class T>
list<T> lorentz_q(list<T> f, T f_pk, T q, T rms);
template<class T, class U>
tuple<list<T>, U, U, U, U, list<T>, list<T>, list<T>, list<T>, T, int, U,
    list<T>, list<T>> calculate_stprod_mono(int nirf_mult, list<T> energy,
    U encomb, U flux_irf, list<T> disk_irf, list<T> gamma_irf, list<T> deltau,
    T min_deltau_frac, int i_rsigmax, list<T> lfreq,  list<T> q, list<T> rms,
    T t_scale);

template <class T>
tuple<T, int> find_nearest(list<T> array, T value) {

    T smallest = numeric_limits<T>::max();
    T small_diff = numeric_limits<T>::max();

    int i = 0;
    int idx = i;

    for (T const &item: array) {
        T diff = abs(item - value);

        if (diff < small_diff) {
            smallest = item;
            small_diff = diff;
            idx = i;
        }

        i++;
    }
    
    return make_tuple(smallest, idx);
}

template <class T>
tuple<list<T>, list<T>, list<T>> calc_dispfrac(list<T> rad, list<T> rad_area,
        T rin, T rcor, T seedff_norm, T seedff_ind, T heatff_norm,
        T heatff_ind) {

    assert (rad.size() == rad_area.size());

    T disptotal = 0;
    list<T> dispfrac(rad.size(), 0);
    list<T> seed_frac_flow(rad.size(), 0);
    list<T> heat_frac_flow(rad.size(), 0);

    typename list<T>::iterator disp, seed, heat, r, area;

    for (disp = dispfrac.begin(), seed = seed_frac_flow.begin(),
            heat = heat_frac_flow.begin(), r = rad.begin(),
            area = rad_area.begin(); disp != dispfrac.end(); disp++, seed++,
            heat++, r++, area++) {
        *disp = *area * pow(*r, -3.0) * (1 - sqrt(rin/(*r)));
        disptotal += *disp;

        if (*r <= rcor) {
            *heat = heatff_norm*pow((*r/rin), heatff_ind);

            if (seedff_norm >= 0) {
                *seed = seedff_norm*pow((*r/rin), heatff_ind);
            } else {
                *seed = 1.0 - *heat;
            }
        }
    }

    for (disp = dispfrac.begin(); disp != dispfrac.end(); disp++) {
        *disp = *disp/disptotal;
    }

    return make_tuple(dispfrac, seed_frac_flow, heat_frac_flow);
}

template <class T, class U>
tuple<list<T>, list<T>, list<T>> calc_illumination_fracs(
    Geometry<list<T>, U> geomod) {

    list<T> omega_cor(geomod.get_omega_cor()),
        frad_disktocor(geomod.get_frad_disktocor()),
        frad_cortodisk (geomod.get_frad_cortodisk_area());

    return make_tuple(omega_cor, frad_disktocor, frad_cortodisk);
}

template <class T>
tuple<list<T>, list<T>, list<T>, list<T>> calc_timing_params(list<T> rad,
        int i_rsigmax, T rcor, int i_rcor, T t_scale, tuple<T, T> disk_tau_par,
        tuple<T, T> cor_tau_par, string lor_model, list<T> lor_par) {
    
    list<T> tau, lfreq, q_list, rms;

    typename list<T>::iterator r, lor = lor_par.begin();
    int i;

    for (r = rad.begin(), i = 0; r != rad.end(); r++, i++) {
        if (*r > rcor) {
            tau.push_back(get<0>(disk_tau_par) * pow((*r/rcor), 
                get<1>(disk_tau_par))*pow(*r, 1.5));
        } else {
            tau.push_back(get<0>(cor_tau_par) * pow((*r/rcor),
                get<1>(cor_tau_par))*pow(*r, 1.5));
        }

        lfreq.push_back(0);
        q_list.push_back(0);
        rms.push_back(0);
    }

    r = rad.begin();

    typename list<T>::iterator t = tau.begin(), freq = lfreq.begin(),
        q = q_list.begin(), rms_value = rms.begin();

    if (lor_model == "continuous") {
        for (i = 0; i < rad.size(); i++) {
            assert(lor_par.size() == 3);
            *next(freq, i) = 1/(*next(t, i)*t_scale);
            *next(q, i) = *next(lor, 0);

            if (i <= i_rsigmax) {
                if (*next(lor, 1) >= 0 && *next(lor, 2) >= 0) {
                    if (*next(r, i) > rcor) {
                        *next(rms_value, i) = sqrt(pow(*next(lor, 1), 2)/
                            (i_rsigmax+1-i_rcor));
                    } else {
                        *next(rms_value, i) =
                            sqrt(pow(*next(lor, 2), 2)/i_rcor);
                    }
                } else if (*next(lor, 1) >= 0 && *next(lor, 2) < 0) {
                    *next(rms_value, i) =
                        sqrt(pow(*next(lor, 1), 2)/(i_rsigmax+1));
                } else if (*next(lor, 1) < 0 && *next(lor, 2) >= 0) {
                    if (*next(r, i) > rcor) {
                        *next(rms_value, i) = -(*next(lor, 1));
                    } else {
                        *next(rms_value, i) =
                            sqrt(pow(*next(lor, 2), 2)/i_rcor);
                    }
                } else {
                    if (*next(r, i) > rcor) {
                        *next(rms_value, i) = -(*next(lor, 1));
                    } else {
                        *next(rms_value, i) = -(*next(lor, 2));
                    }
                }
            }
        }
    }

    if (lor_model == "multi_frequency") {
        for (int j = 0; j < lor_par.size(); j += 3) {
            double ltau = 1/(*next(lor, j)*t_scale);
            if (ltau <= *next(t, tau.size() - 1) && ltau >= *(tau.begin())) {
                auto [ltau2, i] = find_nearest<T>(tau, ltau);
                *next(freq, i) = 1/(*next(t, i) * t_scale);
                *next(q, i) = *next(lor, j+1);
                if (i < i_rsigmax) {
                    *next(rms_value, i) = *next(lor, j+2);
                }
            }
        }
    } else if (lor_model == "multi_radius") {
        for (int j = 0; j < lor_par.size(); j += 3) {
            if (*next(lor, j) <= *next(r, rad.size() - 1) &&
                    *next(lor, j) >= *r) {
                auto [lrad, i] = find_nearest<T>(rad, *next(lor, j));
                *next(freq, i) = 1/(*next(t, i) * t_scale);
                *next(q, i) = *next(lor, j+1);
                if (i <= i_rsigmax) {
                    *next(rms_value, i) = *next(lor, j+2);
                }
            }
        }
    } else if (lor_model == "multi_radius_frequency") {
        for (int j = 0; j < lor_par.size(); j += 4) {
            if (*next(lor, j) <= *next(r, rad.size() - 1) &&
                    *next(lor, j) >= *r) {
                auto [lrad, i] = find_nearest<T>(rad, *next(lor, j));
                if (*next(r, i) > 0) {
                    for (int i_next = 1; i_next < i; i_next++) {
                        if (*next(rms_value, i - i_next) == 0 &&
                                *next(rms_value, i) > 0) {
                            i = i - i_next;
                        }
                    }
                    *next(freq, i) = *next(lor, j+1);
                    *next(q, i) = *next(lor, j+2);
                    if (i < i_rsigmax) {
                        *next(rms_value, i) = *next(lor, j+3);
                    }

                    cout << "Signal radius for " << *next(lor, j+1) <<
                        " Hz is " << *next(r, i) << endl << "with rms " <<
                        *next(rms_value, i) << endl;
                }
            }
        }
    }

    return make_tuple(tau, lfreq, q_list, rms);
}

template <class T>
list<T> calc_propagation_params(list<T> rad, list<T> rad_edge, T rcor,
        tuple<T, T> disk_prop_par, tuple<T, T> cor_prop_par) {
    list<T> deltau;

    typename list<T>::iterator r, edge;

    for (r = rad.begin(), edge = rad_edge.begin(); r != rad.end(); r++,
            edge++) {
        if (*r <= rcor) {
            deltau.push_back(get<0>(cor_prop_par) * (*next(edge, 1) - *edge) * pow(*r, .5) * pow(*r/rcor, get<1>(cor_prop_par)));
        } else {
            deltau.push_back(get<0>(disk_prop_par) * (*next(edge, 1) - *edge) * 
                pow(*r, .5) * pow(*r/rcor, get<1>(disk_prop_par)));
        }
    }

    return deltau;
}

template <class T>
tuple<list<T>, list<T>, list<T>, list<T>, list<T>> calc_radial_time_response(
        list<T> rad, T rcor, list<T> disp_frac, list<T> disktocor_frac,
        list<T> cortodisk_frac, list<T> seed_frac_flow, list<T> heat_frac_flow) {
    
    list<T> ldisk_disp, lseed_disp, lheat, ldisk_rev, lseed_rev;

    typename list<T>::iterator r, disp, fdc, fcd, seed, heat;

    double f_rev = 0, f_return = 0;

    for (fcd = cortodisk_frac.begin(), fdc = disktocor_frac.begin();
            fcd != cortodisk_frac.end(); fcd++, fdc++) {
        f_rev += *fcd; f_return += (*fcd) * (*fdc);
    }

    for (r = rad.begin(), disp = disp_frac.begin(),
            fdc = disktocor_frac.begin(), seed = seed_frac_flow.begin(),
            heat = heat_frac_flow.begin(); r != rad.end(); r++, disp++, fdc++,
            seed++, heat++) {
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
            seed_rev = lseed_rev.begin(); heat != lheat.end(); heat++,
            seed_disp++, seed_rev++) {
        lheat_sum += *heat; lseed_sum += *seed_disp + *seed_rev;
    }

    gamma_mean = get<0>(gamma_par) * pow(lseed_sum/lheat_sum, get<1>(gamma_par));
    u = gamma_mean * get<1>(gamma_par);
    v = gamma_mean - 1;

    for (energy_value = energy.begin(); energy_value != energy.end();
            energy_value++) {
        seed_column.push_back(1-u*log((*energy_value)/e_seed) + (u/v));
        heat_column.push_back(u*log((*energy_value)/e_seed) - (u/v));
    }

    for (seed_disp = lseed_disp.begin(), seed_rev = lseed_rev.begin(),
            heat = lheat.begin(), disk_disp = ldisk_disp.begin(),
            disk_rev = ldisk_rev.begin(); seed_disp != lseed_disp.end();
            seed_disp++, seed_rev++, heat++, disk_disp++, disk_rev++) {
        gamma_irf.push_back(gamma_mean * get<1>(gamma_par) *
            ((*seed_disp + *seed_rev)/lseed_sum - *heat/lheat_sum));
        disk_irf.push_back(*disk_disp + *disk_rev);
        seed_irf.push_back(*seed_disp + *seed_rev);
        seed_div_sum_irf.push_back((*seed_disp + *seed_rev)/lseed_sum);
        heat_div_sum_irf.push_back((*heat)/lheat_sum);
    }

    U firf_seed(seed_column, seed_div_sum_irf),
        firf_heat(heat_column, heat_div_sum_irf);

    U flux_irf = firf_seed + firf_heat;

    return make_tuple(gamma_mean, gamma_irf, flux_irf, disk_irf, seed_irf);
}

template <class T>
tuple<list<T>, T> linear_rebin_irf(T dt, int i_rsigmax, list<T> irf_nbins,
        list<T> irf_binedgefrac, list<T> input_irf, list<T> deltau_scale,
        int nirf) {
    list<T> rebinned_irf(nirf, 0), input_irf2(input_irf.size(), 0);

    typename list<T>::iterator ds = deltau_scale.begin(),
        irf = input_irf.begin(), irf2 = input_irf2.begin(),
        nbins = irf_nbins.begin(), rebinned = rebinned_irf.begin(),
        bef = irf_binedgefrac.begin();
    
    for (int i = i_rsigmax; i >= 0; i--) {
        if (*next(ds, i) > 0) {
            *next(irf2, i) = *next(irf, i) / *next(ds, i);
        }
    }
    
    int irfbin_start;
    int irfbin_stop;

    for (int i = i_rsigmax; i >= 0; i--) {

        T irfbin = 0;

        for (int j = i; j <= i_rsigmax; j++) {
            irfbin += *next(nbins, j);
        }

        if (*next(ds, i) > 0) { 
            irfbin_start = static_cast<int>(irfbin - *next(nbins, i));
            irfbin_stop = static_cast<int>(irfbin_start + *next(nbins, i) - 1);
            for (int j = irfbin_start; j <= irfbin_stop; j++) {
                *next(rebinned, j) = *next(irf2, i);
            }

            if (i < i_rsigmax && *next(bef, i) > 0) {
                *next(rebinned, irfbin_start) = *next(irf2, i) +
                    (*next(bef, i) * (*next(irf2, i+1) - *next(irf2, i)));
            }
            if (i > 0 && *next(bef, i - 1) <= 0) {
                *next(rebinned, irfbin_stop) = *next(irf2, i) +
                    (*next(bef, i-1) * (*next(irf2, i) - *next(irf2, i-1)));
            }
        } else {
            *next(rebinned, irfbin_stop+1) += *next(irf, i);
        }
    }

    *next(rebinned, irfbin_stop+1) /= dt;

    T flux_outer = 0;

    if (i_rsigmax < input_irf.size() - 1) {
        for (int i = i_rsigmax + 1; i < input_irf.size(); i++) {
            flux_outer += *next(irf, i);
        }
    }

    return make_tuple(rebinned_irf, flux_outer);
}

template <class T>
tuple<list<complex<T>>, list<T>, list<T>, list<T>> calc_cross_psd(list<T> freq,
        T dt, list<T> ci_irf, list<T> ref_irf, list <T> irf_nbins,
        list<T> deltau_scale, int i_rsigmax, list<T> f_pk, list<T> q,
        list<T> rms) {
    int nirf = ref_irf.size();

    list<T> wt_pow_spec_ref(nirf/2, 0), wt_pow_spec_ci(nirf/2, 0),
        mod_sig_psd(nirf/2, 0), ref_irf_dum, ci_irf_dum;
    list<complex<T>> wt_cross_spec(nirf/2, 0);

    typename list<T>::iterator ref, ci;

    for (ref = ref_irf.begin(), ci = ci_irf.begin(); ref != ref_irf.end();
            ref++, ci++) {
        ref_irf_dum.push_back((*ref) * dt);
        ci_irf_dum.push_back((*ci) * dt);
    }

    typename list<T>::iterator deltau = deltau_scale.begin(),
        nbins = irf_nbins.begin(), rms_val = rms.begin(),
        f_pk_val = f_pk.begin(), q_val = q.begin(),
        lor, mod_sig, wt_ps_ref, wt_ps_ci;

    typename list<complex<T>>::iterator ref_f, ci_f, wt_cross;

    for (int i = i_rsigmax; i >= 0; i--) {
        if (i == i_rsigmax || *next(deltau, i+1) > 0) {
            if (i < i_rsigmax) {
                T irf = -*next(nbins, i+1);

                for (int j = i+1; j <= i_rsigmax; j++) {
                    irf += *next(nbins, j);
                }

                int irfbin_start = static_cast<int>(irf);
                int irfbin_stop = static_cast<int>(irfbin_start +
                    *next(nbins, i+1) - 1);

                ref = ref_irf_dum.begin(); ci = ci_irf_dum.begin();

                for (int j = irfbin_start; j <= irfbin_stop; j++) {
                    *next(ref, j) = 0;
                    *next(ci, j) = 0;
                }
            }

            if (*next(rms_val, i) > 0) {
                auto lor_psd = lorentz_q<T>(freq, *next(f_pk_val, i),
                    *next(q_val, i), *next(rms_val, i));
                auto ref_fft = FFT<list<T>, T>(ref_irf_dum);
                auto ci_fft = FFT<list<T>, T>(ci_irf_dum);

                for (ref_f = ref_fft.begin(), ci_f = ci_fft.begin(),
                        wt_cross = wt_cross_spec.begin(), lor = lor_psd.begin(),
                        mod_sig = mod_sig_psd.begin(),
                        wt_ps_ref = wt_pow_spec_ref.begin(),
                        wt_ps_ci = wt_pow_spec_ci.begin();
                        ref_f != ref_fft.end(); ref_f++, ci_f++, wt_cross++,
                        lor++, mod_sig++, wt_ps_ref++, wt_ps_ci++) {
                    auto cross = conj(*ci_f) * (*ref_f);
                    *wt_cross += (*lor) * cross;
                    *mod_sig += *lor;

                    *wt_ps_ref += (*lor) * pow(abs(*ref_f),2);
                    *wt_ps_ci += (*lor) * pow(abs(*ci_f),2);
                }
            }
        }
    }

    return make_tuple(wt_cross_spec, wt_pow_spec_ci, wt_pow_spec_ref,
        mod_sig_psd);
}

template <class T, class U>
U lorentz_q(list<T> f, list<T> f_pk, list<T> q, list<T> rms) {
    assert(f_pk.size() == q.size() && q.size() == rms.size());

    U lorentz(f_pk.size(), f.size());

    typename list<T>::iterator f_val, f_pk_val, q_val, rms_val;

    int i, j = 0;

    for (f_pk_val = f_pk.begin(), q_val = q.begin(), rms_val = rms.begin();
            f_pk_val != f_pk.end(); f_pk_val++, q_val++, rms_val++) {
        T f_res = *f_pk_val / sqrt(1+(1/(4*pow(*q_val,2))));
        T r = *rms_val / sqrt(.5-atan(-2*(*q_val))/M_PI);
        T lorentz1 = (1/M_PI)*2*pow(r, 2)*(*q_val)*f_res;

        i = 0;

        for (f_val = f.begin(); f_val != f.end(); f_val++) {
            T lorentz2 = 4*pow(*q_val, 2)*pow(*f_val - f_res, 2);
            lorentz.get_element(j, i) = lorentz1 / (pow(f_res, 2) + lorentz2);

            i++;
        }

        j++;
    }

    return lorentz;
}

template<class T>
list<T> lorentz_q(list<T> f, T f_pk, T q, T rms) {
    T f_res = f_pk / sqrt(1+(1/(4*pow(q, 2))));
    T r = rms / sqrt(.5-atan(-2*q)/M_PI);

    list<T> lorentz(f.size(), 0);

    typename list<T>::iterator lor, f_val;

    for (lor = lorentz.begin(), f_val = f.begin(); lor != lorentz.end(); lor++,
            f_val++) {
        *lor = (1/M_PI*2*pow(r, 2)*q*f_res) / (pow(f_res, 2) +
            (4*pow(q, 2)*pow(*f_val - f_res, 2)));
    }

    return lorentz;
}

template<class T, class U>
tuple<list<T>, U, U, U, U, list<T>, list<T>, list<T>, list<T>, T, int, U,
        list<T>, list<T>> calculate_stprod_mono(int nirf_mult, list<T> energy,
        U encomb, U flux_irf, list<T> disk_irf, list<T> gamma_irf,
        list<T> deltau, T min_deltau_frac, int i_rsigmax, list<T> lfreq,
        list<T> q, list<T> rms, T t_scale) {

    int i = 0;
    T del_min(1E100), deltau_sum_max = 0;
    list<T> deltau_scale;

    typename list<T>::iterator del;

    for (del = deltau.begin(); del != deltau.end(); del++) {
        deltau_scale.push_back(*del * t_scale);
        if (*del > 0 && *del < del_min) {
            del_min = *del;
        }

        if (i < i_rsigmax + 1) { 
            deltau_sum_max += *del * t_scale;
        }

        i++;
    }

    T dt = min_deltau_frac * del_min;
    int nirf = nirf_mult * pow(2, ceil(log2(deltau_sum_max)/dt));

    list<T> irf_nbins(deltau_scale.size(), 0),
        irf_binedgefrac(deltau_scale.size(), 0);
    
    del = deltau_scale.begin();

    int i_irf_max, i_irf_min;
    T del_sum_1 = 0, del_sum_2 = 0, irf_nbins_sum = 0;

    for (i = 0; i <= i_rsigmax; i++) {
        for (int j = 0; j < i; j++) {
            del_sum_1 += *next(del, i);
            del_sum_2 += *next(del, i);
        }
        del_sum_2 += *next(del, i+1);

        i_irf_max = (deltau_sum_max - del_sum_1)/dt - 1;
        i_irf_min = (deltau_sum_max - del_sum_2)/dt;

        irf_nbins.push_back(i_irf_max - i_irf_min + 1);
        irf_nbins_sum += i_irf_max - i_irf_min + 1;
        irf_binedgefrac.push_back((dt*irf_nbins_sum - del_sum_2)/dt);
    }

    U ci_irf(energy.size()+1, nirf);
    list<T> ci_outer(energy.size()+1, 0), ci_mean(energy.size()+1, 1);
    typename list<T>::iterator outer = ci_outer.begin(), mean = ci_mean.begin();

    auto [rebinned, flux_outer] = linear_rebin_irf(dt, i_rsigmax, irf_nbins,
        irf_binedgefrac, disk_irf, deltau_scale, nirf);
    
    *outer = flux_outer;

    typename list<T>::iterator irf = rebinned.begin();
    T ci_sum = 0;

    for (i = 0; i < rebinned.size(); i++) {
        ci_irf.get_element(0, i) = *next(irf, i);
        ci_sum += *next(irf, i);
    }

    *mean = dt * ci_sum + *outer;

    list<T> flux;

    for (i = 1; i <= energy.size(); i++) {
        for (int j = 0; j < get<1>(flux_irf.get_size()); j++) {
            flux.push_back(flux_irf.get_element(i - 1, j));
        }

        auto [rebinned, flux_outer] = linear_rebin_irf(dt, i_rsigmax,
            irf_nbins, irf_binedgefrac, flux, deltau_scale, nirf);

        *next(outer, i) = flux_outer;

        for (int j = 0; j < rebinned.size(); j++) {
            ci_irf.get_element(0, j) = *next(irf, j);
        }
    }

    T minfreq = 1/(dt * nirf);
    list<T> freq;

    for (i = 0; i < nirf/2; i++) {
        freq.push_back(minfreq+i);
    }

    U phlag(get<0>(encomb.get_size()), freq.size()),
        tlag(get<0>(encomb.get_size()), freq.size()),
        psd_ci(get<0>(encomb.get_size()), freq.size()),
        psd_ref(get<0>(encomb.get_size()), freq.size());
    typename list<T>::iterator cross, ps_ci, ps_ref, freq_it = freq.begin();

    list<T> mod_sig_psd;

    for (i = 0; i < get<0>(encomb.get_size()); i++) {
        int j = encomb.get_element(i, 0);
        int k = encomb.get_element(i, 1);

        list<T> irf_j, irf_k;

        for (int l = 0; l < get<1>(ci_irf.get_size()); l++) {
            irf_j.push_back(ci_irf.get_element(j, l));
            irf_k.push_back(ci_irf.get_element(k, l));
        }

        auto [wt_cross_spec, wt_pow_spec_ci, wt_pow_spec_ref, msp] =
            calc_cross_psd(freq, dt, irf_j, irf_k, irf_nbins, deltau_scale,
                i_rsigmax, lfreq, q, rms);

        mod_sig_psd = msp;

        assert(freq.size() == wt_cross_spec.size() == wt_pow_spec_ci.size() == 
            wt_pow_spec_ref.size());

        cross = wt_cross_spec.begin(); ps_ci = wt_pow_spec_ci.begin();
        ps_ref = wt_pow_spec_ref.begin();

        for (int l = 0; l < freq.size(); l++) {
            phlag.get_element(i, l) = arg(*next(cross, l));
            tlag.get_element(i, l) = phlag.get_element(i, l)/
                (2*M_PI*(*next(freq_it, l)));
            psd_ci.get_element(i, l) = *next(ps_ci, l)/pow(*next(mean, j), 2);
            psd_ref.get_element(i, l) = *next(ps_ref, l)/pow(*next(mean, k), 2);
        }
    }

    return make_tuple(freq, phlag, tlag, psd_ci, psd_ref, mod_sig_psd,
        irf_nbins, irf_binedgefrac, deltau_scale, dt, nirf, ci_irf, ci_mean,
        ci_outer);
}

#endif