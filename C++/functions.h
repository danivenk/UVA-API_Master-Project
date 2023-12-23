// make sure only included once
#ifndef functions_H
#define functions_H

// include the standard libraries
#include <vector>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <cassert>
#include <iterator>

// include user libraries
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

// function declarations
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

/**
 * @brief Find the nearest value in an array.
 * 
 * @tparam T is the value type
 * @tparam U is the container type
 * @param array is the array to be searched
 * @param value the value to find the nearest value for
 * @return returns the closest value and it's index
 */
template <typename T, template <typename...> class U>
tuple<T, int> find_nearest(U<T> array, T value) {

    // set the smallest and smallest difference to the max value
    T smallest = numeric_limits<T>::max();
    T small_diff = numeric_limits<T>::max();

    // start with the indices
    int i = 0;
    int idx = i;

    // loop over each item
    for (T const &item: array) {

        // check the difference between the item and value
        T diff = abs(item - value);

        // check if the difference is smaller than the current
        if (diff < small_diff) {

            // save if smaller
            smallest = item;
            small_diff = diff;
            idx = i;
        }

        // increase the index
        i++;
    }
    
    // return the closest value in the array and it's index
    return make_tuple(smallest, idx);
}

/**
 * @brief Calculate the dissipation for each radial bin normalized by the total
 *        dissipation for all radial bins. It also calculates the heat and seed
 *        photon fraction of the dissipated power.
 *        
 *        To allow for unconstrained heating/cooling the `seedff_norm` can be
 *        set as bigger than 0. If it is set to less than 0 to force the total
 *        fraction to be unity.
 *        When `seedff_norm` is set to 0 no seed luminosity is produced inside
 *        the corona.
 * 
 * @tparam T is the value type
 * @tparam U is the container type
 * @param rad radial bins
 * @param rad_area radial area bins
 * @param rin a radius inside the corona
 * @param rcor radius of the corona
 * @param seedff_norm normalization factor for the seed fraction flow
 * @param seedff_ind exponent for the seed fraction flow
 * @param heatff_norm normalization factor for the heat fraction flow
 * @param heatff_ind exponent for the heat fraction flow
 * @return returns the dissipation fraction, seed fraction flow, heat fraction
 *         flow
 */
template <typename T, template <typename...> class U>
tuple<U<T>, U<T>, U<T>> calc_dispfrac(U<T> rad, U<T> rad_area, T rin, T rcor,
        T seedff_norm, T seedff_ind, T heatff_norm, T heatff_ind) {

    // assert the size for the radial bins and the radial area bins as the same
    assert (rad.size() == rad_area.size());

    // set the total dissipation fraction to 0 and make the output arrays with 0
    T disptotal = 0;
    U<T> dispfrac(rad.size(), 0);
    U<T> seed_frac_flow(rad.size(), 0);
    U<T> heat_frac_flow(rad.size(), 0);

    // define the iterators
    typename U<T>::iterator disp, seed, heat, r, area;

    // loop over the arrays to fill the input arrays
    for (disp = dispfrac.begin(), seed = seed_frac_flow.begin(),
            heat = heat_frac_flow.begin(), r = rad.begin(),
            area = rad_area.begin(); disp != dispfrac.end(); disp++, seed++,
            heat++, r++, area++) {

        // calculate the disspation frac for the current bin and add to total
        *disp = *area * pow(*r, -3.0) * (1 - sqrt(rin/(*r)));
        disptotal += *disp;

        // define the heat and seed fraction flows if the bin is outside corona
        if (*r <= rcor) {

            // calculate the heat fraction flow
            *heat = heatff_norm*pow((*r/rin), heatff_ind);

            // set the sum seed and heat to 1 if the norm fractor is negative
            if (seedff_norm >= 0) {
                *seed = seedff_norm*pow((*r/rin), heatff_ind);
            } else {
                *seed = 1.0 - *heat;
            }
        }
    }

    // normalize the disspation fraction
    for (disp = dispfrac.begin(); disp != dispfrac.end(); disp++) {
        *disp = *disp/disptotal;
    }

    // return the disspation fraction, seed fraction flow and heat fraction flow
    return make_tuple(dispfrac, seed_frac_flow, heat_frac_flow);
}

/**
 * @brief Calculates the solid angle & illumination fractions. The illumination
 *        fractions are the fraction of emission that goes from the corona to
 *  	  the disk and the emission that goes from the disk to the corona.
 * 
 * @tparam T is the value type
 * @tparam U is the container type
 * @param geomod is the geometry to calculate the illumination fractions and
 *               solid angle for
 * @return returns the solid angle and illumination fractions for the given
 *         geometry
 */
template <typename T, template <typename...> class U>
tuple<U<T>, U<T>, U<T>> calc_illumination_fracs(Geometry<U, T> geomod) {

    // create the arrays to return from the given geometry
    U<T> omega_cor(geomod.get_omega_cor()),
        frad_disktocor(geomod.get_frad_disktocor()),
        frad_cortodisk (geomod.get_frad_cortodisk_area());

    // return the arrays
    return make_tuple(omega_cor, frad_disktocor, frad_cortodisk);
}

/**
 * @brief Calculates the timing parameters for a given Lorentzian. It these
 *        include the assigned signal power spectrum for each bin from the
 *        Lorentzian, the Lorentzian frequency and the matching Lorentzian
 *        RMS value.
 * 
 *        For this calculation 4 different Lorentzian models exits:
 *        - continuous (Lorentzian signals at all radial bins) ~
 *          Lorentzian model parameters: [q, total_rms_in_disk,
 *          total_rms_in_corona];
 *        - multi_frequency (multiple Lorentzian parameters assigned dependent
 *          on tau) ~
 *          Lorentzian model parameters: [peak_frequency_1, q_1, rms_1,
 *          peak_frequency_2, q_2, rms_2, ...] where len(lor_par) % 3 == 0;
 *        - multi_radius (all the parameters set dependent on the radial bins) ~
 *          Lorentzian model parameters: [radius_1, q_1, rms_1, radius_2, q_2,
 *          rms_2, ...] where len(lor_par) % 3 == 0;
 *        - multi_radius_frequency (is a combination of the two previous) ~
 *          Lorentzian model parameters: [radius_1, peak_frequency_1, q_1,
 *          rms_1, radius_2, peak_frequency_2, q_2, rms_2, ...] where
 *          len(lor_par) % 4 == 0;
 * 
 * @tparam T is the value type
 * @tparam U is the container type
 * @param rad radial bins
 * @param i_rsigmax index of the nearest value for the radial bin for signal max
 * @param rcor radius of the corona
 * @param i_rcor index of the nearest radial bin to the coronal radius
 * @param t_scale scaling factor to convert the propagation delay to seconds
 * @param disk_tau_par parameters for the variability time-scale in the disk 
 * @param cor_tau_par parameters for the variability time-scale in the corona
 * @param lor_model Lorenzian model to use
 * @param lor_par parameters for the chosen Lorentzian
 * @param dbg print out signal radius
 * @return returns the variable time scale, Lorentzian frequency, quality
 *         factor, rms values
 */
template <typename T, template <typename...> class U>
tuple<U<T>, U<T>, U<T>, U<T>> calc_timing_params(U<T> rad,
        int i_rsigmax, T rcor, int i_rcor, T t_scale, tuple<T, T> disk_tau_par,
        tuple<T, T> cor_tau_par, string lor_model, U<T> lor_par, bool dbg) {

    // declare the arrays
    U<T> tau, lfreq, q_list, rms;

    // define the iterators and start the index value
    typename U<T>::iterator r, lor = lor_par.begin();
    int i;

    // loop over the radial bins
    for (r = rad.begin(), i = 0; r != rad.end(); r++, i++) {

        // make sure to only add variabilitt time-scale if outside the corona
        if (*r > rcor) {
            tau.push_back(get<0>(disk_tau_par) * pow((*r/rcor),
                get<1>(disk_tau_par))*pow(*r, 1.5));
        } else {
            tau.push_back(get<0>(cor_tau_par) * pow((*r/rcor),
                get<1>(cor_tau_par))*pow(*r, 1.5));
        }

        // add zeros to the other arrays
        lfreq.push_back(0);
        q_list.push_back(0);
        rms.push_back(0);
    }

    // move the radial bin iterator back to the start
    r = rad.begin();

    // start the iterators for the output arrays
    typename U<T>::iterator t = tau.begin(), freq = lfreq.begin(),
        q = q_list.begin(), rms_value = rms.begin();

    // check for the Lorentz model
    if (lor_model == "continuous") {

        // loop over the indices for the radial bins
        for (i = 0; i < rad.size(); i++) {
            // make sure the amount of Lorentz parameters are a multiple of 3
            assert(lor_par.size() == 3);

            // fill in the frequency and quality factor arrays
            *next(freq, i) = 1/(*next(t, i)*t_scale);
            *next(q, i) = *next(lor, 0);

            // make sure the radial bin is below the maximum signal radius
            if (i <= i_rsigmax) {

                // check the sign of the total rms in the disk and in the corona
                if (*next(lor, 1) >= 0 && *next(lor, 2) >= 0) {
                    // fill the rms values equally over the radial bins
                    if (*next(r, i) > rcor) {
                        *next(rms_value, i) = sqrt(pow(*next(lor, 1), 2)/
                            (i_rsigmax+1-i_rcor));
                    } else {
                        *next(rms_value, i) =
                            sqrt(pow(*next(lor, 2), 2)/i_rcor);
                    }
                } else if (*next(lor, 1) >= 0 && *next(lor, 2) < 0) {
                    // fill the rms values equally over the radial bins
                    *next(rms_value, i) =
                        sqrt(pow(*next(lor, 1), 2)/(i_rsigmax+1));
                } else if (*next(lor, 1) < 0 && *next(lor, 2) >= 0) {
                    // fill the rms values equally over the radial bins for
                    //  inside the corona, for the rest set as the 2nd parameter
                    if (*next(r, i) > rcor) {
                        *next(rms_value, i) = -(*next(lor, 1));
                    } else {
                        *next(rms_value, i) =
                            sqrt(pow(*next(lor, 2), 2)/i_rcor);
                    }
                } else {
                    // fill the the rms values with the values of the parameters
                    //  corresponding to within and outside the corona 
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
        // loop over the sets of 3 parameters
        for (int j = 0; j < lor_par.size(); j += 3) {

            // set the Lorentzian variablility time scale 
            double ltau = 1/(*next(lor, j)*t_scale);

            // check if Lorentzian variability time scale falls with tau value
            if (ltau <= *next(t, tau.size() - 1) && ltau >= *(tau.begin())) {
                // find the nearest Lorenztian value in the tau array
                auto [ltau2, i] = find_nearest<T>(tau, ltau);

                // set the Lorentzian frequency, quality & rms values
                *next(freq, i) = 1/(*next(t, i) * t_scale);
                *next(q, i) = *next(lor, j+1);
                if (i < i_rsigmax) {
                    *next(rms_value, i) = *next(lor, j+2);
                }
            }
        }
    } else if (lor_model == "multi_radius") {
        // loop over the sets of 3 parameters
        for (int j = 0; j < lor_par.size(); j += 3) {

            // check if the Lorenzian parameters are within the radial bin vals
            if (*next(lor, j) <= *next(r, rad.size() - 1) &&
                    *next(lor, j) >= *r) {

                // find the nearest radial bin value to the lorenzian parameter
                auto [lrad, i] = find_nearest<T>(rad, *next(lor, j));
                
                // set the Lorentzian frequency, quality & rms values
                *next(freq, i) = 1/(*next(t, i) * t_scale);
                *next(q, i) = *next(lor, j+1);
                if (i <= i_rsigmax) {
                    *next(rms_value, i) = *next(lor, j+2);
                }
            }
        }
    } else if (lor_model == "multi_radius_frequency") {
        // loop over the sets of 4 parameters
        for (int j = 0; j < lor_par.size(); j += 4) {
            // check if the Lorentzian parameter is within the radial bins
            if (*next(lor, j) <= *next(r, rad.size() - 1) &&
                    *next(lor, j) >= *r) {

                // find the nearest radial bin to the Lorentzian parameter
                auto [lrad, i] = find_nearest<T>(rad, *next(lor, j));

                // make sure the next radial bin is still positive
                if (*next(r, i) > 0) {

                    // loop over the radial bins till the found bin
                    for (int i_next = 1; i_next < i; i_next++) {
                        if (*next(rms_value, i - i_next) == 0 &&
                                *next(rms_value, i) > 0) {
                            i = i - i_next;
                        }
                    }
                    
                    // set the Lorentzian frequency, quality & rms values
                    *next(freq, i) = *next(lor, j+1);
                    *next(q, i) = *next(lor, j+2);
                    if (i < i_rsigmax) {
                        *next(rms_value, i) = *next(lor, j+3);
                    }

                    // print out if dbg bool is set
                    if (dbg) {
                        cout << "Signal radius for " << *next(lor, j+1)
                            << " Hz is " << *next(r, i) << endl << "with rms "
                            << *next(rms_value, i) << endl;
                    }
                }
            }
        }
    }

    // return the variablility time scale, Lorentzian freq, quality & rms arrays
    return make_tuple(tau, lfreq, q_list, rms);
}

/**
 * @brief Calculates the propagation delay across the radial bins according to
 *        the given propagation parameters.
 * 
 * @tparam T is the value type
 * @tparam U is the container type
 * @param rad radial bins
 * @param rad_edge edges of the radial bins
 * @param rcor radius of the corna
 * @param disk_prop_par propagation parameters for the disk
 * @param cor_prop_par propagation parameters for the corona
 * @return returns the propagation delay
 */
template <typename T, template <typename...> class U>
U<T> calc_propagation_params(U<T> rad, U<T> rad_edge, T rcor,
        tuple<T, T> disk_prop_par, tuple<T, T> cor_prop_par) {
    
    // define the array
    U<T> deltau;

    // define the iterators
    typename U<T>::iterator r, edge;

    // loop over the radial bins and their edges
    for (r = rad.begin(), edge = rad_edge.begin(); r != rad.end(); r++,
            edge++) {
        // check if in or outside the corona and calculate the propagation delay
        if (*r <= rcor) {
            deltau.push_back(get<0>(cor_prop_par) * (*next(edge, 1) - *edge) *
                pow(*r, .5) * pow(*r/rcor, get<1>(cor_prop_par)));
        } else {
            deltau.push_back(get<0>(disk_prop_par) * (*next(edge, 1) - *edge) * 
                pow(*r, .5) * pow(*r/rcor, get<1>(disk_prop_par)));
        }
    }

    // return the propagation delay
    return deltau;
}

/**
 * @brief Calculate the disk, seed & heat impulse responses due to the
 *        disspation from each bin.
 * 
 * @tparam T is the value type 
 * @tparam U is the container type
 * @param rad radial bins
 * @param rcor radius of the corona
 * @param disp_frac dissipation fraction
 * @param disktocor_frac luminosity fraction from the disk to the corona
 * @param cortodisk_frac luminosity fraction from the corona to the disk
 * @param thermal_frac fraction of luminosity not turned into thermal energy
 * @param seed_frac_flow fraction of seed photons in the flow
 * @param heat_frac_flow fraction of heat photons in the flow
 * @return returns the luminosities for the disk disspation, seed disspation,
 *  heating luminosity, disk reverberation, seed reverberation
 */
template <typename T, template <typename...> class U>
tuple<U<T>, U<T>, U<T>, U<T>, U<T>> calc_radial_time_response(U<T> rad, T rcor,
        U<T> disp_frac, U<T> disktocor_frac, U<T> cortodisk_frac,
        T thermal_frac, U<T> seed_frac_flow, U<T> heat_frac_flow) {
    
    // define the return arrays
    U<T> ldisk_disp, lseed_disp, lheat, ldisk_rev, lseed_rev;

    // define the iterators
    typename U<T>::iterator r, disp, fdc, fcd, seed, heat;

    // initialize the cor to disk luminosity & cor to disk to cor luminosity
    double f_rev = 0, f_return = 0;

    // add the cor to disk frac and disk to cor frac to the sums
    for (fcd = cortodisk_frac.begin(), fdc = disktocor_frac.begin();
            fcd != cortodisk_frac.end(); fcd++, fdc++) {
        f_rev += *fcd; f_return += (*fcd) * (*fdc);
    }

    // loop over the arrays
    for (r = rad.begin(), disp = disp_frac.begin(),
            fdc = disktocor_frac.begin(), seed = seed_frac_flow.begin(),
            heat = heat_frac_flow.begin(); r != rad.end(); r++, disp++, fdc++,
            seed++, heat++) {

        // add the disk disspation luminosity only if outside the corona
        if (*r > rcor) {
            ldisk_disp.push_back((*disp) * (1 - *fdc));
        } else {
            ldisk_disp.push_back(0);
        }

        // add the seed dissipation luminosity and heating luminosity
        lseed_disp.push_back((*disp) * (*fdc + *seed));
        lheat.push_back((*disp) * (*heat));

        // add the seed and disk reverberation luminosity
        lseed_rev.push_back((thermal_frac) * f_return * ((*disp) *
            (*fdc + *seed) + (*disp) * (*heat))/(1 - f_return));
        ldisk_rev.push_back((thermal_frac) * (f_rev - f_return) *
            ((*disp) * (*fdc + *seed) + (*disp) * (*heat))/(1 - f_return));
    }

    // return the arrays
    return make_tuple(ldisk_disp, lseed_disp, lheat, ldisk_rev, lseed_rev);
}

/**
 * @brief Calculates the impulse response functions (IRFs) for the monochromatic
 *        photon flux, photon index, total disk emission and total seed
 *        emission. For the total emissions they are summed over both the
 *        dissipation and reverberation emissions. The IRFs are calculated
 *        depending on the given luminosities.
 * 
 * @tparam T is the value type
 * @tparam U is the matrix type
 * @tparam V is the container type
 * @param gamma_par parameters for the photon index
 * @param e_seed monochromatic approximation for the seed energy
 * @param energy energy value to use
 * @param ldisk_disp luminosity for the disk dissipation
 * @param lseed_disp luminosity for the seed dissipation
 * @param lheat heating luminosity
 * @param ldisk_rev luminosity for the disk reverberation
 * @param lseed_rev luminosity for the seed reverberation
 * @return returns the photon index mean, impulse response function for the
 *         photon index, photon flux, total disk emission, total seed emission
 */
template <typename T, class U, template <typename...> class V>
tuple<T, V<T>, U, V<T>, V<T>> calc_irfs_mono(tuple<T, T> gamma_par, T e_seed,
        V<T> energy, V<T> ldisk_disp, V<T> lseed_disp, V<T> lheat,
        V<T> ldisk_rev, V<T> lseed_rev) {

    // define the value and array variables
    double lheat_sum, lseed_sum, gamma_mean, u, v;
    V<T> gamma_irf, disk_irf, seed_irf, seed_div_sum_irf, heat_div_sum_irf,
        seed_column, heat_column;

    // define the iterators
    typename V<T>::iterator heat, seed_disp, seed_rev, disk_disp, disk_rev,
        energy_value;

    // set the summation variables to 0
    lheat_sum = 0, lseed_sum = 0;

    // loop over the heating, seed dissipation & seed reverberation luminations
    for (heat = lheat.begin(), seed_disp = lseed_disp.begin(),
            seed_rev = lseed_rev.begin(); heat != lheat.end(); heat++,
            seed_disp++, seed_rev++) {
        // sum the values to the correct sums
        lheat_sum += *heat; lseed_sum += *seed_disp + *seed_rev;
    }

    // find the mean photon indices and the helper variables
    gamma_mean = get<0>(gamma_par) * pow(lseed_sum/lheat_sum, get<1>(gamma_par));
    u = gamma_mean * get<1>(gamma_par);
    v = gamma_mean - 1;

    // loop over the energy array and add them to the columns
    for (energy_value = energy.begin(); energy_value != energy.end();
            energy_value++) {
        seed_column.push_back(1-u*log((*energy_value)/e_seed) + (u/v));
        heat_column.push_back(u*log((*energy_value)/e_seed) - (u/v));
    }

    // loop over the luminosities
    for (seed_disp = lseed_disp.begin(), seed_rev = lseed_rev.begin(),
            heat = lheat.begin(), disk_disp = ldisk_disp.begin(),
            disk_rev = ldisk_rev.begin(); seed_disp != lseed_disp.end();
            seed_disp++, seed_rev++, heat++, disk_disp++, disk_rev++) {
        
        // calculated the impulse response function values
        gamma_irf.push_back(gamma_mean * get<1>(gamma_par) *
            ((*seed_disp + *seed_rev)/lseed_sum - *heat/lheat_sum));
        disk_irf.push_back(*disk_disp + *disk_rev);
        seed_irf.push_back(*seed_disp + *seed_rev);
        seed_div_sum_irf.push_back((*seed_disp + *seed_rev)/lseed_sum);
        heat_div_sum_irf.push_back((*heat)/lheat_sum);
    }

    // create the seed and heat matrices for the flux impulse response function
    U firf_seed(seed_column, seed_div_sum_irf),
        firf_heat(heat_column, heat_div_sum_irf);

    // add the two matrices together
    U flux_irf = firf_seed + firf_heat;

    // return the arrays and matrices
    return make_tuple(gamma_mean, gamma_irf, flux_irf, disk_irf, seed_irf);
}

/**
 * @brief Calculate the luminosities emitted for a given black body energy band
 *        as well as the fraction of photons.
 * 
 * @tparam T is the value type
 * @tparam U is the container type
 * @param disp_frac dissipation fraction
 * @param cortodisk_frac luminosity fraction from the disk to the corona
 * @param disktocor_frac luminosity fraction from the corona to the disk
 * @param lseed_disp luminosity of the seed dissipation
 * @param lheat heating luminosity
 * @param thermal_frac fraction of luminosity not turned into thermal energy
 * @param rcor radius of the corona
 * @param rad radial bins
 * @param rad_area radial area bins
 * @param eband the energy band to use
 * @param kTvar true = non-variable black body; false = variable black body
 * @return returns the band fraction, temperature radial bins, luminosities for this
 *         energy band: disk, seed, disk dissipation, disk reverberation
 */
template <typename T, template <typename...> class U>
tuple<U<T>, U<T>, U<T>, U<T>, U<T>, U<T>> calc_disk_band(U<T> disp_frac,
        U<T> cortodisk_frac, U<T> disktocor_frac, U<T> lseed_disp, U<T> lheat,
        T thermal_frac, T rcor, U<T> rad, U<T> rad_area, tuple<T, T> eband,
        bool kTvar) {

    // define the apery constant ad the black body photon in normalization
    double C_apery = 1.2020569031595942;
    double bb_phintnorm = kTvar ? 1.5 * C_apery * 15/pow(M_PI, 4) :
        2 * C_apery * 15/pow(M_PI, 4);

    // define the iterators
    typename U<T>::iterator disp, fcd, fdc, seed, heat, r, area, kT;

    // define the arrays
    U<T> kT_rad, band_frac, ldiskband_disp, lseedband_disp,
        ldisk_band, lseed_band;

    // define and set the values to 0
    T f_rev = 0, f_return = 0, C_radldisk_rev = 0, radldisk_rev = 0, kT_max = 0,
        kT_val = 0, band = 0, diskband = 0, seedband = 0, diskband_rev_val = 0,
        seedband_rev_val = 0;

    // loop over the arryas
    for (fcd = cortodisk_frac.begin(), fdc = disktocor_frac.begin(),
            heat = lheat.begin(), seed = lseed_disp.begin();
            fcd != cortodisk_frac.end(), fdc != disktocor_frac.end(),
            heat != lheat.end(), seed != lseed_disp.end(); fcd++, fdc++, heat++,
            seed++) {

        // calculate the total fraction from corona to disk
        f_rev += *fcd; f_return += (*fcd) * (*fdc);

        // calculate sum constant for the radial disk reverberation luminosity
        C_radldisk_rev += (*seed + *heat);
    }

    // normalize of the fraction of returning photons
    C_radldisk_rev /= (1-f_return);

    // get the minimum and maximum for the energy band
    auto [E_min, E_max] = eband;

    // creat empty arrays for the disk and seed reverbations
    U<T> ldiskband_rev(rad.size(), 0), lseedband_rev(rad.size(), 0);

    // define he iterators and indices
    typename U<T>::iterator diskband_rev = ldiskband_rev.begin(),
        seedband_rev = lseedband_rev.begin();
    int i = 0, j = 0;

    // loop over the arrays
    for (area = rad_area.begin(), disp = disp_frac.begin(),
            fcd = cortodisk_frac.begin(), r = rad.begin();
            area != rad_area.end(), disp != disp_frac.end(),
            fcd != cortodisk_frac.end(), r != rad.end(); area++, disp++, fcd++,
            r++) {

        // add the values if outside the corona
        if (*r > rcor) {

            // find the radial disk reverberation luminosity
            radldisk_rev = (*fcd)*(thermal_frac)*C_radldisk_rev;

            // find the temperature value and save the max value
            kT_val = pow((*disp + radldisk_rev)/(*area), .25);
            kT_rad.push_back(kT_val);
            kT_max = kT_max < kT_val ? kT_val : kT_max;
        } else {
            kT_val = 0;
            kT_rad.push_back(0);
        }
    }

    // loop over the temperature values and normalize them to the max value
    for (kT = kT_rad.begin(); kT != kT_rad.end(); kT++) { *kT /= kT_max; }

    // loop over the arrays
    for (disp = disp_frac.begin(), fcd = cortodisk_frac.begin(),
            fdc = disktocor_frac.begin(), kT = kT_rad.begin(), r = rad.begin();
            disp != disp_frac.end(), fcd != cortodisk_frac.end(),
            fdc != disktocor_frac.end(), kT != kT_rad.end(), r != rad.end();
            disp++, fcd++, fdc++, kT++, r++) {

        // only add the black body photon flux values outside the corona
        if (*r > rcor) {
            band = bb_phflux<T>(E_min, E_max, *kT, 1000, kTvar) /
                (bb_phintnorm/(*kT));

            band_frac.push_back(band);
        } else {
            band = 0; band_frac.push_back(0);
        }

        // calculate the disk and seed luminosities
        diskband = *disp * (1 - *fdc) * band;
        seedband = *disp * (*fdc) * band;

        // save them to the correct arrays
        ldiskband_disp.push_back(diskband);
        lseedband_disp.push_back(seedband);

        // reset the index
        j = 0;

        // loop over arrays
        for (seed = lseed_disp.begin(), heat = lheat.begin();
                seed != lseed_disp.end(), heat != lheat.end(); seed++, heat++) {
            
            // find the disk and seed reverberation luminosities
            diskband_rev_val = thermal_frac * band * (*fcd * (1 - *fdc)) *
                (*seed + *heat) / (1-f_return);
            seedband_rev_val = thermal_frac * band * (*fcd) * (*fdc) *
                (*seed + *heat) / (1-f_return);

            // save them to the correct arrays
            *next(diskband_rev, j) += diskband_rev_val;
            *next(seedband_rev, j) += seedband_rev_val;

            // increase the index
            j++;
        }

        // increase the index
        i++;
    }

    // make sure the arrays are the same size
    assert(ldiskband_disp.size() == ldiskband_rev.size() &&
        lseedband_disp.size() == lseedband_rev.size());

    // loop over the arrays and sum them to the correct arrays
    for (i = 0; i < ldiskband_disp.size(); i++) {
        ldisk_band.push_back(*next(ldiskband_disp.begin(), i) +
            *next(ldiskband_rev.begin(), i));
        lseed_band.push_back(*next(lseedband_disp.begin(), i) +
            *next(lseedband_rev.begin(), i));
    }

    // return the arrays
    return make_tuple(band_frac, kT_rad, ldisk_band, lseed_band, ldiskband_disp,
        ldiskband_rev);
}

/**
 * @brief Calculates the black body photon flux for the given values.
 * 
 * @tparam T is the value type
 * @param E_min is the minimum energy value to probe at
 * @param E_max is the maximum energy value to probe at
 * @param kT is the temperature
 * @param nE is the number of energy bins
 * @param kTvar true = non-variable black body; false = variable black body
 * @return returns the black body photon flux
 */
template <typename T>
T bb_phflux(T E_min, T E_max, T kT, int nE, bool kTvar) {

    // find the logarithmic stepsize and sets the other values to 0
    T logstep = (log10(E_max) - log10(E_min))/nE, E = 0, dE = 0, bb_flux = 0;

    // create the energy bins and sum the values
    for (int i = 1; i <= nE; i++) {

        // find the energy value and the energy step size
        E = sqrt(pow(10, log10(E_min) + i*logstep) *
            pow(10, log10(E_min) + (i-1)*logstep));
        dE = pow(10, log10(E_min) + i*logstep) -
            pow(10, log10(E_min) + (i-1)*logstep);

        // check for variable black body, calculate and sum
        if (kTvar) {
            bb_flux += nan_to_num(dE*.25*(pow(E, 3)/pow(kT, 5))*exp(E/kT)*
                pow(exp(E/kT) - 1, -2));
        } else {
            bb_flux += nan_to_num(dE*(pow(E, 2)/pow(kT, 4))/(exp(E/kT) - 1));
        }
    }

    // return the flux value
    return bb_flux/(pow(M_PI, 4)/15);
}

/**
 * @brief Linearly rebin the impulse response function (IRF) given. The output
 *        IRF are only filled to the maximum signal radius to save memory.
 * 
 * @tparam T is the value type
 * @tparam U is the container type
 * @param dt is the time step
 * @param i_rsigmax index of the nearest value for the radial bin for signal max
 * @param irf_nbins number of bins for the impulse response function
 * @param irf_binedgefrac fraction at the edge of the bin for the IRF
 * @param input_irf input impulse response function
 * @param deltau_scale scale for the propagation delay
 * @param nirf number of impulse response functions
 * @return returns the rebinned impulse response function & the flux component
 *         outside the maximum singal radius
 */
template <typename T, template <typename...> class U>
tuple<U<T>, T> linear_rebin_irf(T dt, int i_rsigmax, U<T> irf_nbins,
        U<T> irf_binedgefrac, U<T> input_irf, U<T> deltau_scale, int nirf) {
    
    // create empty arrays for the rebinned & converted impulse reponse function
    U<T> rebinned_irf(nirf, 0), input_irf2(input_irf.size(), 0);

    // define the iterators and set them at the start
    typename U<T>::iterator ds = deltau_scale.begin(),
        irf = input_irf.begin(), irf2 = input_irf2.begin(),
        nbins = irf_nbins.begin(), rebinned = rebinned_irf.begin(),
        bef = irf_binedgefrac.begin();

    // convert the input IRF to per time-delay units
    for (int i = i_rsigmax; i >= 0; i--) {
        if (*next(ds, i) > 0) {
            *next(irf2, i) = *next(irf, i) / *next(ds, i);
        }
    }
    
    // set the impulse response fuction bin indices to 0
    int irfbin_start = 0;
    int irfbin_stop = 0;

    // loop from the maximum signal radius index to the start
    for (int i = i_rsigmax; i >= 0; i--) {

        // set the impulse response function bin value to 0
        T irfbin = 0;

        // loop from the maximum signal radius index to the current index
        for (int j = i; j <= i_rsigmax; j++) {
            irfbin += *next(nbins, j);
        }

        // check the sign of the time-delay
        if (*next(ds, i) > 0) {

            // find the index of the start and the end of the rebinned bin
            irfbin_start = static_cast<int>(round(irfbin - *next(nbins, i)
                - 2E-15));
            irfbin_stop = static_cast<int>(round(irfbin_start + *next(nbins, i)
                - 1 - 2E-15));

            // loop from the start to the end of the rebinned bin
            for (int j = irfbin_start; j <= irfbin_stop; j++) {
                // set the value to the converted impulse response function
                *next(rebinned, j) = *next(irf2, i);
            }

            // check the bin edge frac is more than 0 and left of the rsigmax
            if (i < i_rsigmax && *next(bef, i) > 0) {
                *next(rebinned, irfbin_start) = *next(irf2, i) +
                    (*next(bef, i) * (*next(irf2, i+1) - *next(irf2, i)));
            }

            // check index is not 0 and bin edge frac of prev bin is negative
            if (i > 0 && *next(bef, i - 1) <= 0) {
                *next(rebinned, irfbin_stop) = *next(irf2, i) +
                    (*next(bef, i-1) * (*next(irf2, i) - *next(irf2, i-1)));
            }
        } else {
            // sum over the bin next to the stop index
            *next(rebinned, irfbin_stop+1) += *next(irf, i);
        }
    }

    // divide the bin next to the stop index over the time step
    *next(rebinned, irfbin_stop+1) /= dt;

    // set the outer flux to 0
    T flux_outer = 0;

    // check if the bin is right of the maximum signal radius
    if (i_rsigmax < input_irf.size() - 1) {
        // loop over the bins right of maximum signal radius & add to outer flux
        for (int i = i_rsigmax + 1; i < input_irf.size(); i++) {
            flux_outer += *next(irf, i);
        }
    }

    // returns the rebinned irf and the outer flux
    return make_tuple(rebinned_irf, flux_outer);
}

/**
 * @brief Calculate the cross and power spectrum for the given impulse response
 *        function (IRF) and find driving signal power spectrum density (PSD)
 * 
 * @tparam T is the value type 
 * @tparam U is the container type
 * @param freq frequencies
 * @param dt time step
 * @param ci_irf channel of interest impulse response function
 * @param ref_irf reference impulse response function
 * @param irf_nbins number of bins for the impulse response function
 * @param deltau_scale scale for the propagation delay
 * @param i_rsigmax index of the nearest value for the radial bin for signal max
 * @param f_pk peak frequencies
 * @param q quality factor
 * @param rms root mean square
 * @return returns the weighted cross spectrum, signal power spectrum, reference
 *         power spectrum and the signal power spectrum density
 */
template <typename T, template <typename...> class U>
tuple<U<complex<T>>, U<T>, U<T>, U<T>> calc_cross_psd(U<T> freq, T dt,
        U<T> ci_irf, U<T> ref_irf, U<T> irf_nbins, U<T> deltau_scale,
        int i_rsigmax, U<T> f_pk, U<T> q, U<T> rms) {

    // find the size of the impulse response function
    int nirf = ref_irf.size();

    // create the output arrays and 2 dummy arrays
    U<T> wt_pow_spec_ref(nirf/2, 0), wt_pow_spec_ci(nirf/2, 0),
        mod_sig_psd(nirf/2, 0), ref_irf_dum, ci_irf_dum;
    U<complex<T>> wt_cross_spec(nirf/2, 0);

    // define the iterators
    typename U<T>::iterator ref, ci;

    // fill the dummy arrays with the values of the IRFs iterated over timestep
    for (ref = ref_irf.begin(), ci = ci_irf.begin(); ref != ref_irf.end();
            ref++, ci++) {
        ref_irf_dum.push_back((*ref) * dt);
        ci_irf_dum.push_back((*ci) * dt);
    }

    // define some more iterators
    typename U<T>::iterator deltau = deltau_scale.begin(),
        nbins = irf_nbins.begin(), rms_val = rms.begin(),
        f_pk_val = f_pk.begin(), q_val = q.begin(),
        lor, mod_sig, wt_ps_ref, wt_ps_ci;

    // define the iterators for the FFT result and cross spectrum
    typename list<complex<T>>::iterator ref_f, ci_f;
    typename U<complex<T>>::iterator wt_cross;

    // loop from the maximum signal radius to the start
    for (int i = i_rsigmax; i >= 0; i--) {

        // check if currently at the rsigmax or the next deltau is positive
        if (i == i_rsigmax || *next(deltau, i+1) > 0) {

            // make sure under the maximum signal radius
            if (i < i_rsigmax) {

                // find the start and stop for the impulse response function
                T irf = -*next(nbins, i+1);

                for (int j = i+1; j <= i_rsigmax; j++) {
                    irf += *next(nbins, j);
                }

                int irfbin_start = static_cast<int>(round(irf + 2E-16));
                int irfbin_stop = static_cast<int>(round(irfbin_start +
                    *next(nbins, i+1) - 1 + 2E-16));

                // set the iterators to the start
                ref = ref_irf_dum.begin(); ci = ci_irf_dum.begin();

                // set all the values between the start and the stop to 0
                for (int j = irfbin_start; j <= irfbin_stop; j++) {
                    *next(ref, j) = 0;
                    *next(ci, j) = 0;
                }
            }

            // check if the rms value is positive
            if (*next(rms_val, i) > 0) {

                // find the Lorentzian power spectrum density
                auto lor_psd = lorentz_q<T, U>(freq, *next(f_pk_val, i),
                    *next(q_val, i), *next(rms_val, i));

                // find the fourier transform of the dummy IRFs
                auto ref_fft = FFT<U<T>, T>(ref_irf_dum);
                auto ci_fft = FFT<U<T>, T>(ci_irf_dum);

                // loop over the results of the fourier transforms
                for (ref_f = ref_fft.begin(), ci_f = ci_fft.begin(),
                        wt_cross = wt_cross_spec.begin(), lor = lor_psd.begin(),
                        mod_sig = mod_sig_psd.begin(),
                        wt_ps_ref = wt_pow_spec_ref.begin(),
                        wt_ps_ci = wt_pow_spec_ci.begin();
                        ref_f != ref_fft.end(); ref_f++, ci_f++, wt_cross++,
                        lor++, mod_sig++, wt_ps_ref++, wt_ps_ci++) {
                    // find the cross spectrum
                    auto cross = conj(*ci_f) * (*ref_f);

                    // weight the cross spectrum and signal power density
                    *wt_cross += (*lor) * cross;
                    *mod_sig += *lor;

                    // find the weighted power spectra for the signal and ref
                    *wt_ps_ref += (*lor) * pow(abs(*ref_f),2);
                    *wt_ps_ci += (*lor) * pow(abs(*ci_f),2);
                }
            }
        }
    }

    // returns the arrays
    return make_tuple(wt_cross_spec, wt_pow_spec_ci, wt_pow_spec_ref,
        mod_sig_psd);
}

/**
 * @brief Calculate the Lorentzian Signal power density.
 * 
 * @tparam T is the value type
 * @tparam U is the output matrix type
 * @tparam V is the container type
 * @param f frequencies
 * @param f_pk peak frequencies
 * @param q quality value
 * @param rms root mean square
 * @return returns matrix of Lorentzian Signal power densities 
 */
template <typename T, class U, template <typename...> class V>
U lorentz_q(V<T> f, V<T> f_pk, V<T> q, V<T> rms) {
    // make sure the arrays that should be the same size are the same size
    assert(f_pk.size() == q.size() && q.size() == rms.size());

    // create the out matrix
    U lorentz(f_pk.size(), f.size());

    // define the iterators
    typename V<T>::iterator f_val, f_pk_val, q_val, rms_val;

    // set the indices
    int i, j = 0;

    // loop over the values
    for (f_pk_val = f_pk.begin(), q_val = q.begin(), rms_val = rms.begin();
            f_pk_val != f_pk.end(); f_pk_val++, q_val++, rms_val++) {

        // calculate the resonance frequency and normalization constant
        T f_res = *f_pk_val / sqrt(1+(1/(4*pow(*q_val,2))));
        T r = *rms_val / sqrt(.5-atan(-2*(*q_val))/M_PI);

        // find the numerator of the fraction
        T lorentz1 = (1/M_PI)*2*pow(r, 2)*(*q_val)*f_res;

        // set i to 0
        i = 0;

        // loop over the frequencies
        for (f_val = f.begin(); f_val != f.end(); f_val++) {
            // find the Lorentzian function value
            T lorentz2 = 4*pow(*q_val, 2)*pow(*f_val - f_res, 2);
            lorentz.get_element(j, i) = lorentz1 / (pow(f_res, 2) + lorentz2);

            // go to the next i
            i++;
        }

        // go to the next j
        j++;
    }

    // return the found Lorentzian function matrix
    return lorentz;
}

/**
 * @brief Calculate the Lorentzian Signal power density
 * 
 * @tparam T is the value type
 * @tparam U is the container type
 * @param f frequencies
 * @param f_pk peak frequency
 * @param q quality value
 * @param rms root mean square
 * @return returns array of Lorentzian Signal power densities 
 */
template<typename T, template <typename...> class U>
U<T> lorentz_q(U<T> f, T f_pk, T q, T rms) {
    
    // calculate the resonance frequency and normalization constant
    T f_res = f_pk / sqrt(1+(1/(4*pow(q, 2))));
    T r = rms / sqrt(.5-atan(-2*q)/M_PI);

    // create the out array
    U<T> lorentz(f.size(), 0);

    // define the iterators
    typename U<T>::iterator lor, f_val;

    // loop over the frequencies
    for (lor = lorentz.begin(), f_val = f.begin(); lor != lorentz.end(); lor++,
            f_val++) {

        // find the Lorentzian function values
        *lor = (1/M_PI*2*pow(r, 2)*q*f_res) / (pow(f_res, 2) +
            (4*pow(q, 2)*pow(*f_val - f_res, 2)));
    }

    // return the found Lorentzian function array
    return lorentz;
}

/**
 * @brief Calculates the lag frequency, the power spectra for the disk, seed and
 *        the power-law emission from the different bands.
 * 
 * @tparam T is the value type
 * @tparam U is the matrices type 1
 * @tparam V is the matrices type 2
 * @tparam W is the container type
 * @param nirf_mult is a multiplication factor for the IRF bins
 * @param energy energy values to use
 * @param encomb indices for the channel of interest and reference band
 * @param flux_irf impulse response function for the flux
 * @param disk_irf impulse response function for the disk
 * @param gamma_irf impulse response function for the photon index
 * @param deltau propagation delay
 * @param min_deltau_frac minimum fraction for the propagation delay
 * @param i_rsigmax index for the maximum signal radius
 * @param lfreq Lorentzian frequencies
 * @param q quality values
 * @param rms the root mean square values
 * @param t_scale scaling factor to convert the propagation delay to seconds
 * @param dbg print or not print out some values to the console
 * @return returns the frequencies, photon and time lags, power spectrum
 *         densities for the channel of interest and reference IRFs, the signal
 *         power spectrum density, number of bins for the impulse response
 *         function, fraction at the edge of the bin for the impulse response
 *         function, the scaled propagation delays, the time steps, the number
 *         of impulse response function bins, the channel of interest impulse
 *         response function, mean and outer values
 */
template<typename T, class U, class V, template <typename...> class W>
tuple<W<T>, U, U, U, U, W<T>, W<T>, W<T>, W<T>, T, int, U, W<T>, W<T>>
        calculate_stprod_mono(int nirf_mult, W<T> energy, V encomb, U flux_irf,
            W<T> disk_irf, W<T> gamma_irf, W<T> deltau, T min_deltau_frac,
            int i_rsigmax, W<T> lfreq,  W<T> q, W<T> rms, T t_scale, bool dbg) {

    // check if it may print and then print out
    if (dbg) {
        cout << "#######################################" << endl;
        cout << "Calculating mono-energetic spectral-timing products" << endl;
    }

    // set the index, propagation delay minimum and sum
    int i = 0;
    T del_min(1E100), deltau_sum_max = 0;

    // define the scaled propagation delay array
    W<T> deltau_scale;

    // define the iterator
    typename W<T>::iterator del;

    // loop over the propagation delay array, scale it, find min and sum
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

    // define the time step and number of impulse response function bins
    T dt = min_deltau_frac * del_min;
    int nirf = nirf_mult * pow(2, ceil(log2(deltau_sum_max/dt)));

    // check if it may print and then print out
    if (dbg) {
        cout << "Time bin size dt is: " << dt << endl;
        cout << "The maximum propagation delay is: " << deltau_sum_max
            << " and there are " << nirf << " irf bins. " << endl;
    }

    // define the number of IRF bins, edge fraction for IRF bins arrays
    W<T> irf_nbins(deltau_scale.size(), 0),
        irf_binedgefrac(deltau_scale.size(), 0);

    // define some more iterators
    typename W<T>::iterator nbins = irf_nbins.begin(),
        bef = irf_binedgefrac.begin();
    
    // set the iterator at thestart
    del = deltau_scale.begin();

    // define the indices ofthe max and min value of impulse response function
    int i_irf_max, i_irf_min;

    // define the two propagation sums and the IRF number of bins sum
    T del_sum_1, del_sum_2, irf_nbins_sum = 0;

    // loop from the start to the maximum signal radius
    for (i = 0; i <= i_rsigmax; i++) {
        
        // set the sums to 0
        del_sum_1 = 0; del_sum_2 = 0;

        // loop from the start of the propagation delay array to current index 
        for (int j = 0; j < i; j++) {

            // add to the sums
            del_sum_1 += *next(del, j);
            del_sum_2 += *next(del, j);
        }
        
        // include the current index to the second sum
        del_sum_2 += *next(del, i);

        // set the indices for the impulse response function min and max
        i_irf_max =
            static_cast<int>(round((deltau_sum_max - del_sum_1)/dt)) - 1;
        i_irf_min =
            static_cast<int>(round((deltau_sum_max - del_sum_2)/dt));

        // set the values to the arrays and sum the number if IRF bins value
        *next(nbins, i) = i_irf_max - i_irf_min + 1;
        irf_nbins_sum += i_irf_max - i_irf_min + 1;
        *next(bef, i) = (dt*irf_nbins_sum - del_sum_2)/dt;
    }

    // create the channel of interest IRF matrix, outer & mean array
    U ci_irf(energy.size()+1, nirf);
    W<T> ci_outer(energy.size()+1, 0), ci_mean(energy.size()+1, 1);

    // define the iterators for those arrays and set them
    typename W<T>::iterator outer = ci_outer.begin(), mean = ci_mean.begin();

    // rebind the disk impulse response function
    auto [rebinned, flux_outer] = linear_rebin_irf<T, W>(dt, i_rsigmax,
        irf_nbins, irf_binedgefrac, disk_irf, deltau_scale, nirf);
    
    // set the first outer value to the flux from the rebinning
    *outer = flux_outer;

    // define the iterator for the rebinned irf and the channel of interest sum
    typename W<T>::iterator irf = rebinned.begin();
    T ci_sum = 0;

    // loop over the rebinned irf array and add the value to the matrix and sum
    for (i = 0; i < rebinned.size(); i++) {
        ci_irf.get_element(0, i) = *next(irf, i);
        ci_sum += *next(irf, i);
    }

    // set the first mean vaue
    *mean = dt * ci_sum + *outer;

    // define the flux array and it's iterator
    W<T> flux(get<1>(flux_irf.get_size()), 0);
    typename W<T>::iterator f = flux.begin();

    // loop over the energies
    for (i = 1; i <= energy.size(); i++) {
        // loop over the flux IRF values and add the IRF value to the array
        for (int j = 0; j < get<1>(flux_irf.get_size()); j++) {
            *next(f, j) = flux_irf.get_element(i - 1, j);
        }

        // rebin the flux values
        auto [rebinned, flux_outer] = linear_rebin_irf<T, W>(dt, i_rsigmax,
            irf_nbins, irf_binedgefrac, flux, deltau_scale, nirf);

        // reset the iterator to the begin of the array
        irf = rebinned.begin();

        // set the first outer value to the flux from the rebinning
        *next(outer, i) = flux_outer;

        // loop over the rebinned impulse response function and set the ci_irf
        for (int j = 0; j < rebinned.size(); j++) {
            ci_irf.get_element(i, j) = *next(irf, j);
        }
    }

    // define the minimum frequency and defie the frequency array
    T minfreq = 1/(dt * nirf);
    W<T> freq;

    // set the frequency values to the array for half of the number of irf bins
    for (i = 1; i <= nirf/2; i++) {
        freq.push_back(minfreq*i);
    }

    // define photon lag, time lag, psd channel of interest & reference matrices
    U phlag(get<0>(encomb.get_size()), freq.size()),
        tlag(get<0>(encomb.get_size()), freq.size()),
        psd_ci(get<0>(encomb.get_size()), freq.size()),
        psd_ref(get<0>(encomb.get_size()), freq.size());
    
    // define some more iterators
    typename W<T>::iterator ps_ci, ps_ref, freq_it = freq.begin();
    typename W<complex<T>>::iterator cross;

    // define the signal power spectrum density array
    W<T> mod_sig_psd;

    // loop over the indices for the channel of interest and reference band
    for (i = 0; i < get<0>(encomb.get_size()); i++) {

        // get the indices from the encomb matrix
        int j = encomb.get_element(i, 0);
        int k = encomb.get_element(i, 1);

        // check if it may print and then print out
        if (dbg) {
            cout << "CI mean, ref mean, CI outer, ref outer : "
                << *next(mean, j) << " " << *next(mean, k) << " "
                << *next(outer, j) << " " << " " << *next(outer, k) << endl;
        }

        // create the impulse response functions for each index
        W<T> irf_j, irf_k;

        // loop over ci_IRF matrix and add to their respective IRF arrays
        for (int l = 0; l < get<1>(ci_irf.get_size()); l++) {
            irf_j.push_back(ci_irf.get_element(j, l));
            irf_k.push_back(ci_irf.get_element(k, l));
        }

        // calculate the cross & power spectra
        auto [wt_cross_spec, wt_pow_spec_ci, wt_pow_spec_ref, msp] =
            calc_cross_psd<T, W>(freq, dt, irf_j, irf_k, irf_nbins,
                deltau_scale, i_rsigmax, lfreq, q, rms);

        // save the power spectrum densitys
        mod_sig_psd = msp;

        // make sure the arrays are the same size
        assert(freq.size() == wt_cross_spec.size() && wt_pow_spec_ci.size() == 
            wt_pow_spec_ref.size() && wt_cross_spec.size() ==
            wt_pow_spec_ci.size());

        // save the cross and power spectra
        cross = wt_cross_spec.begin(); ps_ci = wt_pow_spec_ci.begin();
        ps_ref = wt_pow_spec_ref.begin();

        // loop over the frequencies
        for (int l = 0; l < freq.size(); l++) {
            // calculate the photon lag and time lag, add them to their matrices
            phlag.get_element(i, l) = arg(*next(cross, l));
            tlag.get_element(i, l) = phlag.get_element(i, l)/
                (2*M_PI*(*next(freq_it, l)));

            // calculate power spectrum densities and add them to their matrices
            psd_ci.get_element(i, l) = *next(ps_ci, l)/pow(*next(mean, j), 2);
            psd_ref.get_element(i, l) = *next(ps_ref, l)/pow(*next(mean, k), 2);
        }

        // check if it may print and then print out
        if (dbg) {
            cout << "Calculated for energies [";
            
            // get the size of the encomb along the first axis
            auto lmax = get<1>(encomb.get_size());

            // loop for this many items
            for (int l = 0; l < lmax; l++) {

                // print out all the values
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

    // return the arrays, matrices and values
    return make_tuple(freq, phlag, tlag, psd_ci, psd_ref, mod_sig_psd,
        irf_nbins, irf_binedgefrac, deltau_scale, dt, nirf, ci_irf, ci_mean,
        ci_outer);
}

/**
 * @brief Return max value for infinity and 0 for NaN. Returns value if neither
 * 
 * @tparam T is the value type
 * @param val is the value to be checked
 * @return is the converted value
 */
template<typename T>
T nan_to_num(T val) {
    return isinf(val) ? numeric_limits<T>::max() : isnan(val) ? 0 : val;
}

#endif