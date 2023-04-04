// include the geometry base class
#include "geometry.h"

// make sure it is only included once
#ifndef an_sphere_H
#define an_shpere_H

/**
 * @brief Describes an analytical sphere geometry that is based on the base
 *        geometry and can also be used in the calc_illumination_fracs function.
 *        It calculates the solid angle, corona to disk fraction & disk to
 *        corona fraction for each radial bin.
 * 
 * @tparam T is the container type
 * @tparam U is the value type
 */
template <template <typename...> class T, typename U>
class AN_Sphere : public Geometry<T, U> {
    public:
        /**
         * @brief Constructor for an AN_Sphere class
         * 
         * @param r radial bins;
         * @param r_area radial area bins;
         * @param r_cor radius of the corona;
         */
        AN_Sphere(T<U> r, T<U> r_area, U r_cor) :
                Geometry<T, U>(r, r_area, r_cor) {
            setup_geometry();
            this->initialize_area_lists();
        };

    private:

        /**
         * @brief Set up the geometry
         */
        void setup_geometry() {
            
            // define the iterators
            typename T<U>::iterator fcd, fdc, r;

            // loop over the arrays to be filled and the radial bins
            for (fcd = this->_frad_cortodisk.begin(),
                    fdc = this->_frad_disktocor.begin(), r = this->_r.begin();
                    fcd != this->_frad_cortodisk.end(),
                    fdc != this->_frad_disktocor.end(), r != this->_r.end();
                    fcd++, fdc++, r++) {
                
                // only add if the radius is bigger than the corona radius
                if (*r > this->_r_cor) {
                    
                    // calculate the disk to corona and corona to disk fraction
                    *fdc = (asin(this->_r_cor/(*r)) - 
                        sqrt(pow((this->_r_cor/(*r)), 2) -
                        pow((this->_r_cor/(*r)), 4)))/M_PI;
                    *fcd = *fdc/(2*M_PI*pow(this->_r_cor, 2));
                }
            }
        }
};

#endif