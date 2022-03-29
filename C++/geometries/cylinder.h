#include "geometry.h"

template <class T, class U>
class Cylinder : public Geometry<T, U> {
    public:
        Cylinder(T r, T r_area, U parms) : Geometry<T, U>(r, r_area, parms) {
            setup_geometry();
            this->initialize_area_lists();
        };

    private:
        void setup_geometry() {
            
            typename T::iterator omega, fcd, fdc, r;

            assert(tuple_size<U>{} == 4);

            auto [r_cor, h_cor, nz, nphi] = this->_parms;

            double dz = h_cor/nz, dphi = 2*M_PI/nphi;

            T x_arr, y_arr, z_arr, da_arr, ypos_arr, zpos_arr;

            for (double z = 0; z < h_cor; z += dz) {
                for (double phi = 0; phi < 2*M_PI; phi += dphi) {
                    double x = cos(phi)*r_cor;
                    double y = sin(phi)*r_cor;
                    double da = dz*dphi*r_cor;

                    x_arr.push_back(x); y_arr.push_back(y); z_arr.push_back(z);
                    da_arr.push_back(x);

                    double ypos = -y;
                    double zpos = -z;
                    
                    ypos_arr.push_back(ypos); zpos_arr.push_back(zpos);
                }
            }

            typename T::iterator x, y, z, da, ypos, zpos;

            for (omega = this->_omega_cor.begin(),
                    fcd = this->_frad_cortodisk.begin(),
                    fdc = this->_frad_disktocor.begin(), r = this->_r.begin();
                    omega != this->_omega_cor.end(),
                    fcd != this->_frad_cortodisk.end(),
                    fdc != this->_frad_disktocor.end(), r != this->_r.end();
                    omega++, fcd++, fdc++, r++) {
                double omegawt = 0;

                for (x = x_arr.begin(), y = y_arr.begin(), z = z_arr.begin(),
                        da = da_arr.begin(), ypos = ypos_arr.begin(),
                        zpos = zpos_arr.begin(); x != x_arr.end(),
                        y != y_arr.end(), z != z_arr.end(), da != da_arr.end(),
                        ypos != ypos_arr.end(), zpos != zpos_arr.end(); x++,
                        y++, z++, da++, ypos++, zpos++) {
                    double xpos = (*r)-(*x);
                    double rpos = pow(pow(xpos, 2) + pow(*ypos,2) +
                                        pow(*zpos, 2), .5);

                    xpos /= rpos;
                    *ypos /= rpos;
                    *zpos /= rpos;

                    double daomega = (*da)*((*x)*(xpos) + (*y)*(*ypos))/
                        (pow(rpos, 2)*r_cor);
                    double daomegawt = daomega*abs(*zpos);

                    if (daomega > 0 && *r > r_cor) {
                        *omega += daomega;
                    }
                    if (daomegawt > 0 && *r > r_cor) {
                        omegawt += daomegawt;
                    }
                }

                *fdc = omegawt/M_PI;
                *fcd = *fdc/((2*M_PI*r_cor*h_cor) + (M_PI*pow(r_cor, 2)));
            }
        }
};