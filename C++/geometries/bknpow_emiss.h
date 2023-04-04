// include the geometry base class
#include "geometry.h"

// make sure it is only included once
#ifndef bknpow_emiss_H
#define bknpow_emiss_H

/**
 * @brief Describes a broken powerlaw emission geometry that is based on the
 *        base geometry and can also be used in the calc_illumination_fracs
 *        function.
 *        It calculates the solid angle, corona to disk fraction & disk to
 *        corona fraction for each radial bin.
 * 
 * @tparam T is the container type
 * @tparam U is the value type
 */
template <template <typename...> class T, typename U>
class BKNpow_Emiss : public Geometry<T, U> {
    public:
        /**
         * @brief Constructor for a BKNpow_Emiss class
         * 
         * @param r radial bins;
         * @param r_area radial area bins;
         * @param r_cor radius of the corona;
         * @param r_bk1 radius of break point 1;
         * @param r_bk2 radius of break point 2;
         * @param em_ind1 induced emission of the first section;
         * @param em_ind2 induced emission of the second section;
         * @param em_ind3 induced emission of the thrid section;
         * @param total_cortodisk total corona to disk fraction;
         */
        BKNpow_Emiss(T<U> r, T<U> r_area, U r_cor, U r_bk1, U r_bk2, U em_ind1,
                U em_ind2, U em_ind3, U total_cortodisk) :
                Geometry<T, U>(r, r_area, r_cor) {
            _r_bk1 = r_bk1, _r_bk2 = r_bk2, _em_ind1 = em_ind1,
            _em_ind2 = em_ind2, _em_ind3 = em_ind3, _tot_ctd = total_cortodisk;
            setup_geometry();
            this->initialize_area_lists();
            normalize_geometry();
        };

    private:
        // first set the total fraction to 0
        double total_frac = 0;
        // define the vaiables for the values
        U _r_bk1, _r_bk2, _em_ind1, _em_ind2, _em_ind3, _tot_ctd;

        /**
         * @brief Set up the geometry
         */
        void setup_geometry() {
            
            // define the iterators
            typename T<U>::iterator omega, fcd, fdc, r, area;

            // loop over the arrays to be filled and the radial bins
            for (omega = this->_omega_cor.begin(),
                    fcd = this->_frad_cortodisk.begin(),
                    fdc = this->_frad_disktocor.begin(), r = this->_r.begin(),
                    area = this->_r_area.begin();
                    omega != this->_omega_cor.end(),
                    fcd != this->_frad_cortodisk.end(),
                    fdc != this->_frad_disktocor.end(), r != this->_r.end(),
                    area != this->_r_area.end(); omega++, fcd++, fdc++, r++,
                    area++) {
                // check the section the radius falls in make sure it is inside
                //  corona and calculate the corona to disk fraction
                if (*r < _r_bk1 && *r > this->_r_cor) {
                    *fcd = pow((*r/_r_bk1), _em_ind1);
                } else if (*r >= _r_bk1 && *r < _r_bk2 && *r > this->_r_cor) {
                    *fcd = pow((*r/_r_bk1), _em_ind2);
                } else if (*r >= _r_bk2 && *r > this->_r_cor) {
                    *fcd = pow((_r_bk2/_r_bk1), _em_ind2) * 
                        pow((*r/_r_bk2), _em_ind3);
                }

                // calculate the disk to corona fraction
                *fdc = *fcd*2*M_PI*pow(this->_r_cor, 2);

                // sum the total fraction
                total_frac += (*fcd)*(*area);
            }
        }

        /**
         * @brief Normalize the geometry
         */
        void normalize_geometry() {

            // define the iterators
            typename T<U>::iterator fcd, fdc;

            // loop over the arrays to be normalized
            for (fcd = this->_frad_cortodisk_area.begin(),
                    fdc = this->_frad_disktocor.begin();
                    fcd != this->_frad_cortodisk_area.end(),
                    fdc != this->_frad_disktocor.end();
                    fcd++, fdc++) {
                
                // normalize both arrays
                *fcd = (*fcd) * _tot_ctd / total_frac;
                *fdc = (*fdc) * _tot_ctd / total_frac;
            }
        }
};

#endif