// include the geometry base class
#include "geometry.h"

// make sure it is only included once
#ifndef sphere_H
#define sphere_H

/**
 * @brief Describe a sphere geometry that is based on the base geometry and can
 *        also be used in the calc_illumination_fracs function.
 *        It calculates the solid angle, corona to disk fraction & disk to
 *        corona fraction for each radial bin.
 * 
 * @tparam T is the container type
 * @tparam U is the value type
 */
template <template <typename...> class T, typename U>
class Sphere : public Geometry<T, U> {
    public:
        /**
         * @brief Constructor for a Sphere class
         * 
         * @param r radial bins;
         * @param r_area radial area bins;
         * @param r_cor radius of the corona;
         * @param nphi number of grid elements in the phi direction;
         * @param ntheta number of grid elements in the theta direction;
         */
        Sphere(T<U> r, T<U> r_area, U r_cor, int nphi, int ntheta) :
                Geometry<T, U>(r, r_area, r_cor) {
            _nphi = nphi, _ntheta = ntheta;
            setup_geometry();
            this->initialize_area_lists();
        };

    private:
        // number of gird elements in the theta & phi direction;
        int _ntheta, _nphi;

        /**
         * @brief Set up the geometry
         */
        void setup_geometry() {
            
            // define the iterators
            typename T<U>::iterator omega, fcd, fdc, r;

            // set the steps in the theta & phi directions
            double dtheta = M_PI/(2*_ntheta), dphi = 2*M_PI/_nphi;

            // define the arrays
            T<U> x_arr, y_arr, z_arr, da_arr, ypos_arr, zpos_arr;

            // define the current value of theta & phi
            double theta_val = dtheta/2, phi_val = dphi/2;

            // loop over each grid value
            for (int i = 0; i < _ntheta; i++) {
                for (int j = 0; j < _nphi; j++) {

                    // define the x, y, z values and the step in area
                    double x = sin(theta_val)*cos(phi_val)*this->_r_cor;
                    double y = sin(theta_val)*sin(phi_val)*this->_r_cor;
                    double z = cos(theta_val)*this->_r_cor;
                    double da = pow(this->_r_cor, 2)*dtheta*dphi*sin(theta_val);

                    // add the values to the correct arrays
                    x_arr.push_back(x); y_arr.push_back(y); z_arr.push_back(z);
                    da_arr.push_back(da);

                    // define the y & z positions
                    double ypos = -y;
                    double zpos = -z;
                    
                    // add the positions to the correct arrays
                    ypos_arr.push_back(ypos); zpos_arr.push_back(zpos);

                    // increase the phi value
                    phi_val += dphi;
                }
                // increase the theta value
                theta_val += dtheta;
            }

            // define some more iterators for the next loop
            typename T<U>::iterator x, y, z, da, ypos, zpos;

            // loop over the arrays to be filled and the radial bins
            for (omega = this->_omega_cor.begin(),
                    fcd = this->_frad_cortodisk.begin(),
                    fdc = this->_frad_disktocor.begin(), r = this->_r.begin();
                    omega != this->_omega_cor.end(),
                    fcd != this->_frad_cortodisk.end(),
                    fdc != this->_frad_disktocor.end(), r != this->_r.end();
                    omega++, fcd++, fdc++, r++) {
                // first set the weighted solid angle to 0
                double omegawt = 0;

                // then loop over the values calculated earlier
                for (x = x_arr.begin(), y = y_arr.begin(), z = z_arr.begin(),
                        da = da_arr.begin(), ypos = ypos_arr.begin(),
                        zpos = zpos_arr.begin(); x != x_arr.end(),
                        y != y_arr.end(), z != z_arr.end(), da != da_arr.end(),
                        ypos != ypos_arr.end(), zpos != zpos_arr.end(); x++,
                        y++, z++, da++, ypos++, zpos++) {
                    
                    // calculate the x position & the radius position
                    double xpos = (*r)-(*x);
                    double rpos = pow(pow(xpos, 2) + pow(*ypos, 2) +
                                      pow(*zpos, 2), .5);

                    // calculate the solid angle and the weighted solid angle
                    double daomega = (*da)*((*x)*(xpos/rpos) +
                        (*y)*(*ypos/rpos) + (*z)*(*zpos/rpos))/
                        (pow(rpos, 2)*this->_r_cor);
                    double daomegawt = daomega*abs(*zpos/rpos);

                    // only sum if within the radius of the corona and positive
                    if (daomega > 0 && *r > this->_r_cor) {
                        *omega += daomega;
                    }
                    if (daomegawt > 0 && *r > this->_r_cor) {
                        omegawt += daomegawt;
                    }
                }

                // calculate the disk to corona and corona to disk fraction
                *fdc = omegawt/M_PI;
                *fcd = *fdc/(2*M_PI*pow(this->_r_cor, 2));
            }
        }
};

#endif