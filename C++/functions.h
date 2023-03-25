#ifndef functions_H
#define functions_H

#include <vector>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <cassert>
#include <iterator>

#include "includes.h"
#include "geometries/geometry.h"
#include "geometries/bknpow_emiss.h"
#include "geometries/an_sphere.h"
#include "geometries/sphere.h"
#include "geometries/cylinder.h"
#include "geometries/inv_cone.h"
#include "geometries/piecewise_emiss.h"
#include "FFT/FFT.h"
#include "array.h"

template <typename T, template <typename...> class U>
tuple<T, int> find_nearest(U<T> array, T value);
template <typename T, template <typename...> class U>
tuple<U<T>, U<T>, U<T>> calc_dispfrac(U<T> rad, U<T> rad_area, T rin, T rcor,
    T seedff_norm, T seedff_ind, T heatff_norm, T heatff_ind);
template <typename T, template <typename...> class U>
tuple<U<T>, U<T>, U<T>> calc_illumination_fracs(Geometry<U, T> geomod);
template <typename T, template <typename...> class U>
tuple<U<T>, U<T>, U<T>, U<T>> calc_timing_params(U<T> rad, int i_rsigmax,
    T rcor, int i_rcor, T t_scale, tuple<T, T> disk_tau_par,
    tuple<T, T> cor_tau_par, string lor_model, U<T> lor_par, bool dbg=false);
template <typename T, template <typename...> class U>
U<T> calc_propagation_parms(U<T> rad, U<T> rad_edge, T rcor,
    tuple<T, T> disk_prop_par, tuple<T, T> cor_prop_par);
template <typename T, template <typename...> class U>
tuple<U<T>, U<T>, U<T>, U<T>, U<T>> calc_radial_time_response(U<T> rad, T rcor,
    U<T> disp_frac, U<T> disktocor_frac, U<T> cortodisk_frac, T thermal_frac,
    U<T> seed_frac_flow, U<T> heat_frac_flow);
template <typename T, class U, template <typename...> class V>
tuple<T, V<T>, U, V<T>, V<T>> calc_irfs_mono(tuple<T, T> gamma_par, T e_seed,
    V<T> energy, V<T> ldisk_disp, V<T> lseed_disp, V<T> lheat, V<T> ldisk_rev,
    V<T> lseed_rev);
template <typename T, template <typename...> class U>
tuple<U<T>, U<T>, U<T>, U<T>, U<T>, U<T>> calc_disk_band(U<T> disp_frac,
    U<T> cortodisk_frac, U<T> disktocor_frac, U<T> lseed_disp, U<T> lheat,
    T thermal_frac, T rcor, U<T> rad, U<T> rad_area, tuple<T, T> eband,
    bool kTvar);
template <typename T>
T bb_phflux(T E_min, T E_max, T kT, int nE, bool kTvar);
template <typename T, template <typename...> class U>
tuple<U<T>, T> linear_rebin_irf(T dt, int i_rsigmax, U<T> irf_nbins,
    U<T> irf_binedgefrac, U<T> input_irf, U<T> deltau_scale, int nirf);
template <typename T, template <typename...> class U>
tuple<U<complex<T>>, U<T>, U<T>, U<T>> calc_cross_psd(U<T> freq, T dt,
    U<T> ci_irf, U<T> ref_irf, U<T> irf_nbins, U<T> deltau_scale,int i_rsigmax,
    U<T> f_pk, U<T> q, U<T> rms);
template <typename T, class U, template <typename...> class V>
U lorentz_q(V<T> f, V<T> f_pk, V<T> q, V<T> rms);
template<class T, template <typename...> class U>
U<T> lorentz_q(U<T> f, T f_pk, T q, T rms);
template<typename T, class U, class V, template <typename...> class W>
tuple<W<T>, U, U, U, U, W<T>, W<T>, W<T>, W<T>, T, int, U, W<T>, W<T>>
    calculate_stprod_mono(int nirf_mult, W<T> energy, V encomb, U flux_irf,
        W<T> disk_irf, W<T> gamma_irf, W<T> deltau, T min_deltau_frac,
        int i_rsigmax, W<T> lfreq,  W<T> q, W<T> rms, T t_scale, bool dbg = false);
template<typename T>
T nan_to_num(T val);

template <typename T, template <typename...> class U>
tuple<T, int> find_nearest(U<T> array, T value) {

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

template <typename T, template <typename...> class U>
tuple<U<T>, U<T>, U<T>> calc_dispfrac(U<T> rad, U<T> rad_area, T rin, T rcor,
        T seedff_norm, T seedff_ind, T heatff_norm, T heatff_ind) {

    assert (rad.size() == rad_area.size());

    T disptotal = 0;
    U<T> dispfrac(rad.size(), 0);
    U<T> seed_frac_flow(rad.size(), 0);
    U<T> heat_frac_flow(rad.size(), 0);

    typename U<T>::iterator disp, seed, heat, r, area;

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

template <typename T, template <typename...> class U>
tuple<U<T>, U<T>, U<T>> calc_illumination_fracs(Geometry<U, T> geomod) {

    U<T> omega_cor(geomod.get_omega_cor()),
        frad_disktocor(geomod.get_frad_disktocor()),
        frad_cortodisk (geomod.get_frad_cortodisk_area());

    return make_tuple(omega_cor, frad_disktocor, frad_cortodisk);
}

template <typename T, template <typename...> class U>
tuple<U<T>, U<T>, U<T>, U<T>> calc_timing_params(U<T> rad,
        int i_rsigmax, T rcor, int i_rcor, T t_scale, tuple<T, T> disk_tau_par,
        tuple<T, T> cor_tau_par, string lor_model, U<T> lor_par, bool dbg) {

    U<T> tau, lfreq, q_list, rms;

    typename U<T>::iterator r, lor = lor_par.begin();
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

    typename U<T>::iterator t = tau.begin(), freq = lfreq.begin(),
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

                    if (dbg) {
                        cout << "Signal radius for " << *next(lor, j+1)
                            << " Hz is " << *next(r, i) << endl << "with rms "
                            << *next(rms_value, i) << endl;
                    }
                }
            }
        }
    }

    return make_tuple(tau, lfreq, q_list, rms);
}

template <typename T, template <typename...> class U>
U<T> calc_propagation_params(U<T> rad, U<T> rad_edge, T rcor,
        tuple<T, T> disk_prop_par, tuple<T, T> cor_prop_par) {
    U<T> deltau;

    typename U<T>::iterator r, edge;

    for (r = rad.begin(), edge = rad_edge.begin(); r != rad.end(); r++,
            edge++) {
        if (*r <= rcor) {
            deltau.push_back(get<0>(cor_prop_par) * (*next(edge, 1) - *edge) *
                pow(*r, .5) * pow(*r/rcor, get<1>(cor_prop_par)));
        } else {
            deltau.push_back(get<0>(disk_prop_par) * (*next(edge, 1) - *edge) * 
                pow(*r, .5) * pow(*r/rcor, get<1>(disk_prop_par)));
        }
    }

    return deltau;
}

template <typename T, template <typename...> class U>
tuple<U<T>, U<T>, U<T>, U<T>, U<T>> calc_radial_time_response(U<T> rad, T rcor,
        U<T> disp_frac, U<T> disktocor_frac, U<T> cortodisk_frac,
        T thermal_frac, U<T> seed_frac_flow, U<T> heat_frac_flow) {
    
    U<T> ldisk_disp, lseed_disp, lheat, ldisk_rev, lseed_rev;

    typename U<T>::iterator r, disp, fdc, fcd, seed, heat;

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

        lseed_disp.push_back((*disp) * (*fdc + *seed));
        lheat.push_back((*disp) * (*heat));

        lseed_rev.push_back((thermal_frac) * f_return * ((*disp) *
            (*fdc + *seed) + (*disp) * (*heat))/(1 - f_return));
        ldisk_rev.push_back((thermal_frac) * (f_rev - f_return) *
            ((*disp) * (*fdc + *seed) + (*disp) * (*heat))/(1 - f_return));
    }

    return make_tuple(ldisk_disp, lseed_disp, lheat, ldisk_rev, lseed_rev);
}

template <typename T, class U, template <typename...> class V>
tuple<T, V<T>, U, V<T>, V<T>> calc_irfs_mono(tuple<T, T> gamma_par, T e_seed,
        V<T> energy, V<T> ldisk_disp, V<T> lseed_disp, V<T> lheat,
        V<T> ldisk_rev, V<T> lseed_rev) {

    double lheat_sum, lseed_sum, gamma_mean, u, v;
    V<T> gamma_irf, disk_irf, seed_irf, seed_div_sum_irf, heat_div_sum_irf,
        seed_column, heat_column;

    typename V<T>::iterator heat, seed_disp, seed_rev, disk_disp, disk_rev,
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

template <typename T, template <typename...> class U>
tuple<U<T>, U<T>, U<T>, U<T>, U<T>, U<T>> calc_disk_band(U<T> disp_frac,
        U<T> cortodisk_frac, U<T> disktocor_frac, U<T> lseed_disp, U<T> lheat,
        T thermal_frac, T rcor, U<T> rad, U<T> rad_area, tuple<T, T> eband,
        bool kTvar) {
    double C_apery = 1.2020569031595942;
    double bb_phintnorm = kTvar ? 1.5 * C_apery * 15/pow(M_PI, 4) :
        2 * C_apery * 15/pow(M_PI, 4);

    typename U<T>::iterator disp, fcd, fdc, seed, heat, r, area, kT;
    U<T> kT_rad, band_frac, ldiskband_disp, lseedband_disp,
        ldisk_band, lseed_band;
    T f_rev = 0, f_return = 0, C_radldisk_rev = 0, radldisk_rev = 0, kT_max = 0,
        kT_val = 0, band = 0, diskband = 0, seedband = 0, diskband_rev_val = 0,
        seedband_rev_val = 0;

    for (fcd = cortodisk_frac.begin(), fdc = disktocor_frac.begin(),
            heat = lheat.begin(), seed = lseed_disp.begin();
            fcd != cortodisk_frac.end(), fdc != disktocor_frac.end(),
            heat != lheat.end(), seed != lseed_disp.end(); fcd++, fdc++, heat++,
            seed++) {
        f_rev += *fcd; f_return += (*fcd) * (*fdc);
        C_radldisk_rev += (*seed + *heat);
    }

    C_radldisk_rev /= (1-f_return);
    auto [E_min, E_max] = eband;

    U<T> ldiskband_rev(rad.size(), 0), lseedband_rev(rad.size(), 0);
    typename U<T>::iterator diskband_rev = ldiskband_rev.begin(),
        seedband_rev = lseedband_rev.begin();
    int i = 0, j = 0;

    for (area = rad_area.begin(), disp = disp_frac.begin(),
            fcd = cortodisk_frac.begin(), r = rad.begin();
            area != rad_area.end(), disp != disp_frac.end(),
            fcd != cortodisk_frac.end(), r != rad.end(); area++, disp++, fcd++,
            r++) {
        if (*r > rcor) {
            radldisk_rev = (*fcd)*(thermal_frac)*C_radldisk_rev;
            kT_val = pow((*disp + radldisk_rev)/(*area), .25);
            kT_rad.push_back(kT_val);
            kT_max = kT_max < kT_val ? kT_val : kT_max;
        } else {
            kT_val = 0;
            kT_rad.push_back(0);
        }
    }

    for (kT = kT_rad.begin(); kT != kT_rad.end(); kT++) { *kT /= kT_max; }

    for (disp = disp_frac.begin(), fcd = cortodisk_frac.begin(),
            fdc = disktocor_frac.begin(), kT = kT_rad.begin(), r = rad.begin();
            disp != disp_frac.end(), fcd != cortodisk_frac.end(),
            fdc != disktocor_frac.end(), kT != kT_rad.end(), r != rad.end();
            disp++, fcd++, fdc++, kT++, r++) {
        if (*r > rcor) {
            band = bb_phflux<T>(E_min, E_max, *kT, 1000, kTvar) /
                (bb_phintnorm/(*kT));

            band_frac.push_back(band);
        } else {
            band = 0; band_frac.push_back(0);
        }

        diskband = *disp * (1 - *fdc) * band;
        seedband = *disp * (*fdc) * band;

        ldiskband_disp.push_back(diskband);
        lseedband_disp.push_back(seedband);

        j = 0;

        for (seed = lseed_disp.begin(), heat = lheat.begin();
                seed != lseed_disp.end(), heat != lheat.end(); seed++, heat++) {
            diskband_rev_val = thermal_frac * band * (*fcd * (1 - *fdc)) *
                (*seed + *heat) / (1-f_return);
            seedband_rev_val = thermal_frac * band * (*fcd) * (*fdc) *
                (*seed + *heat) / (1-f_return);

            *next(diskband_rev, j) += diskband_rev_val;
            *next(seedband_rev, j) += seedband_rev_val;
            j++;
        }

        i++;
    }

    assert(ldiskband_disp.size() == ldiskband_rev.size() &&
        lseedband_disp.size() == lseedband_rev.size());

    for (i = 0; i < ldiskband_disp.size(); i++) {
        ldisk_band.push_back(*next(ldiskband_disp.begin(), i) +
            *next(ldiskband_rev.begin(), i));
        lseed_band.push_back(*next(lseedband_disp.begin(), i) +
            *next(lseedband_rev.begin(), i));
    }

    return make_tuple(band_frac, kT_rad, ldisk_band, lseed_band, ldiskband_disp,
        ldiskband_rev);
}

template <typename T>
T bb_phflux(T E_min, T E_max, T kT, int nE, bool kTvar) {
    T logstep = (log10(E_max) - log10(E_min))/nE, E = 0, dE = 0, bb_flux = 0;

    for (int i = 1; i <= nE; i++) {
        E = sqrt(pow(10, log10(E_min) + i*logstep) *
            pow(10, log10(E_min) + (i-1)*logstep));
        dE = pow(10, log10(E_min) + i*logstep) -
            pow(10, log10(E_min) + (i-1)*logstep);

        if (kTvar) {
            bb_flux += nan_to_num(dE*.25*(pow(E, 3)/pow(kT, 5))*exp(E/kT)*
                pow(exp(E/kT) - 1, -2));
        } else {
            bb_flux += nan_to_num(dE*(pow(E, 2)/pow(kT, 4))/(exp(E/kT) - 1));
        }
    }

    return bb_flux/(pow(M_PI, 4)/15);
}

template <typename T, template <typename...> class U>
tuple<U<T>, T> linear_rebin_irf(T dt, int i_rsigmax, U<T> irf_nbins,
        U<T> irf_binedgefrac, U<T> input_irf, U<T> deltau_scale, int nirf) {
    U<T> rebinned_irf(nirf, 0), input_irf2(input_irf.size(), 0);

    typename U<T>::iterator ds = deltau_scale.begin(),
        irf = input_irf.begin(), irf2 = input_irf2.begin(),
        nbins = irf_nbins.begin(), rebinned = rebinned_irf.begin(),
        bef = irf_binedgefrac.begin();

    for (int i = i_rsigmax; i >= 0; i--) {
        if (*next(ds, i) > 0) {
            *next(irf2, i) = *next(irf, i) / *next(ds, i);
        }
    }
    
    int irfbin_start = 0;
    int irfbin_stop = 0;

    for (int i = i_rsigmax; i >= 0; i--) {

        T irfbin = 0;

        for (int j = i; j <= i_rsigmax; j++) {
            irfbin += *next(nbins, j);
        }

        if (*next(ds, i) > 0) { 
            irfbin_start = static_cast<int>(round(irfbin - *next(nbins, i)
                - 2E-15));
            irfbin_stop = static_cast<int>(round(irfbin_start + *next(nbins, i)
                - 1 - 2E-15));

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

template <typename T, template <typename...> class U>
tuple<U<complex<T>>, U<T>, U<T>, U<T>> calc_cross_psd(U<T> freq, T dt,
        U<T> ci_irf, U<T> ref_irf, U<T> irf_nbins, U<T> deltau_scale,
        int i_rsigmax, U<T> f_pk, U<T> q, U<T> rms) {
    int nirf = ref_irf.size();

    U<T> wt_pow_spec_ref(nirf/2, 0), wt_pow_spec_ci(nirf/2, 0),
        mod_sig_psd(nirf/2, 0), ref_irf_dum, ci_irf_dum;
    U<complex<T>> wt_cross_spec(nirf/2, 0);

    typename U<T>::iterator ref, ci;

    for (ref = ref_irf.begin(), ci = ci_irf.begin(); ref != ref_irf.end();
            ref++, ci++) {
        ref_irf_dum.push_back((*ref) * dt);
        ci_irf_dum.push_back((*ci) * dt);
    }

    typename U<T>::iterator deltau = deltau_scale.begin(),
        nbins = irf_nbins.begin(), rms_val = rms.begin(),
        f_pk_val = f_pk.begin(), q_val = q.begin(),
        lor, mod_sig, wt_ps_ref, wt_ps_ci;

    typename list<complex<T>>::iterator ref_f, ci_f;
    typename U<complex<T>>::iterator wt_cross;

    for (int i = i_rsigmax; i >= 0; i--) {
        if (i == i_rsigmax || *next(deltau, i+1) > 0) {
            if (i < i_rsigmax) {
                T irf = -*next(nbins, i+1);

                for (int j = i+1; j <= i_rsigmax; j++) {
                    irf += *next(nbins, j);
                }

                int irfbin_start = static_cast<int>(round(irf + 2E-16));
                int irfbin_stop = static_cast<int>(round(irfbin_start +
                    *next(nbins, i+1) - 1 + 2E-16));

                ref = ref_irf_dum.begin(); ci = ci_irf_dum.begin();

                for (int j = irfbin_start; j <= irfbin_stop; j++) {
                    *next(ref, j) = 0;
                    *next(ci, j) = 0;
                }
            }

            if (*next(rms_val, i) > 0) {
                auto lor_psd = lorentz_q<T, U>(freq, *next(f_pk_val, i),
                    *next(q_val, i), *next(rms_val, i));
                auto ref_fft = FFT<U<T>, T>(ref_irf_dum);
                auto ci_fft = FFT<U<T>, T>(ci_irf_dum);

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

template <typename T, class U, template <typename...> class V>
U lorentz_q(V<T> f, V<T> f_pk, V<T> q, V<T> rms) {
    assert(f_pk.size() == q.size() && q.size() == rms.size());

    U lorentz(f_pk.size(), f.size());

    typename V<T>::iterator f_val, f_pk_val, q_val, rms_val;

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

template<typename T, template <typename...> class U>
U<T> lorentz_q(U<T> f, T f_pk, T q, T rms) {
    T f_res = f_pk / sqrt(1+(1/(4*pow(q, 2))));
    T r = rms / sqrt(.5-atan(-2*q)/M_PI);

    U<T> lorentz(f.size(), 0);

    typename U<T>::iterator lor, f_val;

    for (lor = lorentz.begin(), f_val = f.begin(); lor != lorentz.end(); lor++,
            f_val++) {
        *lor = (1/M_PI*2*pow(r, 2)*q*f_res) / (pow(f_res, 2) +
            (4*pow(q, 2)*pow(*f_val - f_res, 2)));
    }

    return lorentz;
}

template<typename T, class U, class V, template <typename...> class W>
tuple<W<T>, U, U, U, U, W<T>, W<T>, W<T>, W<T>, T, int, U, W<T>, W<T>>
        calculate_stprod_mono(int nirf_mult, W<T> energy, V encomb, U flux_irf,
            W<T> disk_irf, W<T> gamma_irf, W<T> deltau, T min_deltau_frac,
            int i_rsigmax, W<T> lfreq,  W<T> q, W<T> rms, T t_scale, bool dbg) {

    if (dbg) {
        cout << "#######################################" << endl;
        cout << "Calculating mono-energetic spectral-timing products" << endl;
    }

    int i = 0;
    T del_min(1E100), deltau_sum_max = 0;
    W<T> deltau_scale;

    typename W<T>::iterator del;

    for (del = deltau.begin(); del != deltau.end(); del++) {
        deltau_scale.push_back(*del * t_scale);
        if (*del * t_scale > 0 && *del * t_scale < del_min) {
            del_min = *del * t_scale;
        }

        if (i < i_rsigmax + 1) { 
            deltau_sum_max += *del * t_scale;
        }

        i++;
    }

    T dt = min_deltau_frac * del_min;
    int nirf = nirf_mult * pow(2, ceil(log2(deltau_sum_max/dt)));

    if (dbg) {
        cout << "Time bin size dt is: " << dt << endl;
        cout << "The maximum propagation delay is: " << deltau_sum_max
            << " and there are " << nirf << " irf bins. " << endl;
    }

    W<T> irf_nbins(deltau_scale.size(), 0),
        irf_binedgefrac(deltau_scale.size(), 0);

    typename W<T>::iterator nbins = irf_nbins.begin(),
        bef = irf_binedgefrac.begin();
    
    del = deltau_scale.begin();

    int i_irf_max, i_irf_min;
    T del_sum_1, del_sum_2, irf_nbins_sum = 0;

    for (i = 0; i <= i_rsigmax; i++) {
        del_sum_1 = 0; del_sum_2 = 0;
        for (int j = 0; j < i; j++) {
            del_sum_1 += *next(del, j);
            del_sum_2 += *next(del, j);
        }
        del_sum_2 += *next(del, i);

        i_irf_max =
            static_cast<int>(round((deltau_sum_max - del_sum_1)/dt)) - 1;
        i_irf_min =
            static_cast<int>(round((deltau_sum_max - del_sum_2)/dt));

        *next(nbins, i) = i_irf_max - i_irf_min + 1;
        irf_nbins_sum += i_irf_max - i_irf_min + 1;
        *next(bef, i) = (dt*irf_nbins_sum - del_sum_2)/dt;
    }

    U ci_irf(energy.size()+1, nirf);
    W<T> ci_outer(energy.size()+1, 0), ci_mean(energy.size()+1, 1);
    typename W<T>::iterator outer = ci_outer.begin(), mean = ci_mean.begin();

    auto [rebinned, flux_outer] = linear_rebin_irf<T, W>(dt, i_rsigmax,
        irf_nbins, irf_binedgefrac, disk_irf, deltau_scale, nirf);
    
    *outer = flux_outer;

    typename W<T>::iterator irf = rebinned.begin();
    T ci_sum = 0;

    for (i = 0; i < rebinned.size(); i++) {
        ci_irf.get_element(0, i) = *next(irf, i);
        ci_sum += *next(irf, i);
    }

    *mean = dt * ci_sum + *outer;

    W<T> flux(get<1>(flux_irf.get_size()), 0);
    typename W<T>::iterator f = flux.begin();

    for (i = 1; i <= energy.size(); i++) {
        for (int j = 0; j < get<1>(flux_irf.get_size()); j++) {
            *next(f, j) = flux_irf.get_element(i - 1, j);
        }

        auto [rebinned, flux_outer] = linear_rebin_irf<T, W>(dt, i_rsigmax,
            irf_nbins, irf_binedgefrac, flux, deltau_scale, nirf);

        irf = rebinned.begin();

        *next(outer, i) = flux_outer;

        for (int j = 0; j < rebinned.size(); j++) {
            ci_irf.get_element(i, j) = *next(irf, j);
        }
    }

    T minfreq = 1/(dt * nirf);
    W<T> freq;

    for (i = 1; i <= nirf/2; i++) {
        freq.push_back(minfreq*i);
    }

    U phlag(get<0>(encomb.get_size()), freq.size()),
        tlag(get<0>(encomb.get_size()), freq.size()),
        psd_ci(get<0>(encomb.get_size()), freq.size()),
        psd_ref(get<0>(encomb.get_size()), freq.size());
    typename W<T>::iterator ps_ci, ps_ref, freq_it = freq.begin();
    typename W<complex<T>>::iterator cross;

    W<T> mod_sig_psd;

    for (i = 0; i < get<0>(encomb.get_size()); i++) {
        int j = encomb.get_element(i, 0);
        int k = encomb.get_element(i, 1);

        if (dbg) {
            cout << "CI mean, ref mean, CI outer, ref outer : "
                << *next(mean, j) << " " << *next(mean, k) << " "
                << *next(outer, j) << " " << " " << *next(outer, k) << endl;
        }

        W<T> irf_j, irf_k;

        for (int l = 0; l < get<1>(ci_irf.get_size()); l++) {
            irf_j.push_back(ci_irf.get_element(j, l));
            irf_k.push_back(ci_irf.get_element(k, l));
        }

        auto [wt_cross_spec, wt_pow_spec_ci, wt_pow_spec_ref, msp] =
            calc_cross_psd<T, W>(freq, dt, irf_j, irf_k, irf_nbins,
                deltau_scale, i_rsigmax, lfreq, q, rms);

        mod_sig_psd = msp;

        assert(freq.size() == wt_cross_spec.size() && wt_pow_spec_ci.size() == 
            wt_pow_spec_ref.size() && wt_cross_spec.size() ==
            wt_pow_spec_ci.size());

        cross = wt_cross_spec.begin(); ps_ci = wt_pow_spec_ci.begin();
        ps_ref = wt_pow_spec_ref.begin();

        for (int l = 0; l < freq.size(); l++) {
            phlag.get_element(i, l) = arg(*next(cross, l));
            tlag.get_element(i, l) = phlag.get_element(i, l)/
                (2*M_PI*(*next(freq_it, l)));
            psd_ci.get_element(i, l) = *next(ps_ci, l)/pow(*next(mean, j), 2);
            psd_ref.get_element(i, l) = *next(ps_ref, l)/pow(*next(mean, k), 2);
        }

        if (dbg) {
            cout << "Calculated for energies [";
            
            auto lmax = get<1>(encomb.get_size());

            for (int l = 0; l < lmax; l++) {
                cout << encomb.get_element(i, l);
                if (l != lmax - 1) {
                    cout << " ";
                }
                else {
                    cout << "]" << endl;
                }
            }
        }
    }

    return make_tuple(freq, phlag, tlag, psd_ci, psd_ref, mod_sig_psd,
        irf_nbins, irf_binedgefrac, deltau_scale, dt, nirf, ci_irf, ci_mean,
        ci_outer);
}

template<typename T>
T nan_to_num(T val) {
    return isinf(val) ? numeric_limits<T>::max() : isnan(val) ? 0 : val;
}

#endif