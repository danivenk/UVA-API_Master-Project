#include "geometry.h"

template <template <typename...> class T, typename U>
class Cylinder : public Geometry<T, U> {
    public:
        Cylinder(T<U> r, T<U> r_area, U r_cor, U h_cor, int nz, int nphi) :
                Geometry<T, U>(r, r_area, r_cor) {
            _h_cor = h_cor, _nz = nz, _nphi = nphi;
            setup_geometry();
            this->initialize_area_lists();
        };

    private:
        U _h_cor;
        int _nz, _nphi;

        void setup_geometry() {
            
            typename T<U>::iterator omega, fcd, fdc, r;

            double dz = _h_cor/_nz, dphi = 2*M_PI/_nphi;

            T<U> x_arr, y_arr, z_arr, da_arr, ypos_arr, zpos_arr;

            for (double z = 0; z/dz < _nz; z += dz) {
                for (double phi = 0; phi/dphi < _nphi; phi += dphi) {
                    double x = cos(phi)*this->_r_cor;
                    double y = sin(phi)*this->_r_cor;
                    double da = dz*dphi*this->_r_cor;

                    x_arr.push_back(x); y_arr.push_back(y); z_arr.push_back(z);
                    da_arr.push_back(da);

                    double ypos = -y;
                    double zpos = -z;
                    
                    ypos_arr.push_back(ypos); zpos_arr.push_back(zpos);
                }
            }

            typename T<U>::iterator x, y, z, da, ypos, zpos;

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

                    double daomega = (*da)*((*x)*(xpos/rpos) +
                        (*y)*(*ypos/rpos))/(pow(rpos, 2)*this->_r_cor);
                    double daomegawt = daomega*abs(*zpos/rpos);

                    if (daomega > 0 && *r > this->_r_cor) {
                        *omega += daomega;
                    }
                    if (daomegawt > 0 && *r > this->_r_cor) {
                        omegawt += daomegawt;
                    }
                }

                *fdc = omegawt/M_PI;
                *fcd = *fdc/((2*M_PI*this->_r_cor*_h_cor) +
                    (M_PI*pow(this->_r_cor, 2)));
            }
        }
};