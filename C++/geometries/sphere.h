#include "geometry.h"

template <template <typename...> class T, typename U>
class Sphere : public Geometry<T, U> {
    public:
        Sphere(T<U> r, T<U> r_area, U r_cor, int nphi, int ntheta) :
                Geometry<T, U>(r, r_area, r_cor) {
            _nphi = nphi, _ntheta = ntheta;
            setup_geometry();
            this->initialize_area_lists();
        };

    private:
        int _ntheta, _nphi;

        void setup_geometry() {
            
            typename T<U>::iterator omega, fcd, fdc, r;

            double dtheta = M_PI/(2*_ntheta), dphi = 2*M_PI/_nphi;

            T<U> x_arr, y_arr, z_arr, da_arr, ypos_arr, zpos_arr;

            double theta_val = dtheta/2, phi_val = dphi/2;

            for (int i = 0; i < _ntheta; i++) {
                for (int j = 0; j < _nphi; j++) {

                    double x = sin(theta_val)*cos(phi_val)*this->_r_cor;
                    double y = sin(theta_val)*sin(phi_val)*this->_r_cor;
                    double z = cos(theta_val)*this->_r_cor;
                    double da = pow(this->_r_cor, 2)*dtheta*dphi*sin(theta_val);

                    x_arr.push_back(x); y_arr.push_back(y); z_arr.push_back(z);
                    da_arr.push_back(da);

                    double ypos = -y;
                    double zpos = -z;
                    
                    ypos_arr.push_back(ypos); zpos_arr.push_back(zpos);

                    phi_val += dphi;
                }
                theta_val += dtheta;
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
                    double rpos = pow(pow(xpos, 2) + pow(*ypos, 2) +
                                      pow(*zpos, 2), .5);

                    double daomega = (*da)*((*x)*(xpos/rpos) +
                        (*y)*(*ypos/rpos) + (*z)*(*zpos/rpos))/
                        (pow(rpos, 2)*this->_r_cor);
                    double daomegawt = daomega*abs(*zpos/rpos);

                    if (daomega > 0 && *r > this->_r_cor) {
                        *omega += daomega;
                    }
                    if (daomegawt > 0 && *r > this->_r_cor) {
                        omegawt += daomegawt;
                    }
                }

                *fdc = omegawt/M_PI;
                *fcd = *fdc/(2*M_PI*pow(this->_r_cor, 2));
            }
        }
};
