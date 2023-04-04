// include the geometry base class
#include "geometry.h"

// make sure it is only included once
#ifndef inv_cone_H
#define inv_cone_H

/**
 * @brief Describes a inverted cone geometry that is based on the base geometry
 *        and can also be used in the calc_illumination_fracs function.
 *        It calculates the solid angle, corona to disk fraction & disk to
 *        corona fraction for each radial bin.
 * 
 * @tparam T is the container type
 * @tparam U is the value type
 */
template <template <typename...> class T, typename U>
class Inv_Cone : public Geometry<T, U> {
    public:
        /**
         * @brief Constructor for a Inv_Cone class
         * 
         * @param r radial bins;
         * @param r_area radial area bins;
         * @param r_cor radius of the corona;
         * @param h_cor height of the corona;
         * @param r_top radius at the top;
         * @param nphi number of grid elements in the phi direction;
         * @param nz number of grid elements in the z direction;
         */
        Inv_Cone(T<U> r, T<U> r_area, U r_cor, U h_cor, U r_top, int nphi,
                int nz) : Geometry<T, U>(r, r_area, r_cor) {
            _h_cor = h_cor, _r_top = r_top, _nphi = nphi, _nz = nz;
            setup_geometry();
            this->initialize_area_lists();
        };

    private:
        // height of the corona, radius at the top;
        U _h_cor, _r_top;
        // number of grid elements in the phi & z direction;
        int _nz, _nphi;

        /**
         * @brief Set up the geometry
         */
        void setup_geometry() {
            
            // define the iterators
            typename T<U>::iterator omega, fcd, fdc, r;

            // set the steps in the z & phi directions
            double dz = _h_cor/_nz, dphi = 2*M_PI/_nphi;

            // define the arrays
            T<U> r_cone_arr, da_arr, x_arr, y_arr, z_arr,
                ypos_arr, zpos_arr;
            
            // define the cone angle and current value of z & phi
            double cone_angle = atan((_r_top-this->_r_cor)/_h_cor);
            double z_val = dz/2, phi_val = dphi/2;

            // loop over each grid value
            for (int i = 0; i < _nz; i++) {
                for (int j = 0; j < _nphi; j++) {

                    // define the radius of the cone and the step in area
                    double r_cone = this->_r_cor +
                        z_val*(_r_top-this->_r_cor)/_h_cor;
                    double da = dphi*r_cone*dz/cos(cone_angle);

                    // define the x, y, z values of the current element
                    double x = cos(phi_val)*r_cone;
                    double y = sin(phi_val)*r_cone;
                    double _z = -r_cone*(_r_top-this->_r_cor)/_h_cor;

                    // add the values to the correct arrays
                    r_cone_arr.push_back(r_cone); da_arr.push_back(da);
                    x_arr.push_back(x); y_arr.push_back(y); z_arr.push_back(_z);

                    // define the y & z positions
                    double ypos = -y;
                    double zpos = -z_val;
                    
                    // add the positions to the correct arrays
                    ypos_arr.push_back(ypos); zpos_arr.push_back(zpos);

                    // increase the phi value
                    phi_val += dphi;
                }
                // increase the z value
                z_val += dz;
            }

            // define some more iterators for the next for loop
            typename T<U>::iterator r_cone, da, x, y, z, ypos, zpos;

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
                for (r_cone = r_cone_arr.begin(), da = da_arr.begin(),
                        x = x_arr.begin(), y = y_arr.begin(), z = z_arr.begin(),
                        ypos = ypos_arr.begin(), zpos = zpos_arr.begin();
                        r_cone != r_cone_arr.end(), da != da_arr.end(),
                        x != x_arr.end(), y != y_arr.end(), z != z_arr.end(),
                        ypos != ypos_arr.end(), zpos != zpos_arr.end();
                        r_cone++, da++, x++, y++, z++, ypos++, zpos++) {

                    // calculate the x position & the radius position
                    double xpos = (*r)-(*x);
                    double rpos = pow(pow(xpos, 2) + pow(*ypos,2) +
                                        pow(*zpos, 2), .5);

                    // calculate the radial position of the corona
                    double rcorpos = pow(pow(*r_cone, 2) + pow(*z, 2), .5);

                    // calculate the solid angle and the weighted solid angle
                    double daomega = (*da)*((*x)*(xpos/rpos) +
                        (*y)*(*ypos/rpos) + (*z)*(*zpos/rpos))/
                        (pow(rpos, 2)*rcorpos);
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
                *fcd = *fdc/(M_PI*(pow(_r_top, 2) + (this->_r_cor + _r_top) *
                    pow(pow((_r_top - this->_r_cor),2)+pow(_h_cor, 2), .5)));
            }
        }
};

#endif