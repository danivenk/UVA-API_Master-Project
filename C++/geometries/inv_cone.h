#include "geometry.h"

template <template <typename...> class T, typename U>
class Inv_Cone : public Geometry<T, U> {
    public:
        Inv_Cone(T<U> r, T<U> r_area, U r_cor, U h_cor, U r_top, int nphi,
                int nz) : Geometry<T, U>(r, r_area, r_cor) {
            _h_cor = h_cor, _r_top = r_top, _nphi = nphi, _nz = nz;
            setup_geometry();
            this->initialize_area_lists();
        };

    private:
        U _h_cor, _r_top;
        int _nz, _nphi;

        void setup_geometry() {
            
            typename T<U>::iterator omega, fcd, fdc, r;

            double dz = _h_cor/_nz, dphi = 2*M_PI/_nphi;

            T<U> r_cone_arr, da_arr, x_arr, y_arr, z_arr,
                ypos_arr, zpos_arr;
                
            double cone_angle = atan((_r_top-this->_r_cor)/_h_cor);
            double z_val = dz/2, phi_val = dphi/2;

            for (int i = 0; i < _nz; i++) {
                for (int j = 0; j < _nphi; j++) {
                    double r_cone = this->_r_cor +
                        z_val*(_r_top-this->_r_cor)/_h_cor;
                    double da = dphi*r_cone*dz/cos(cone_angle);

                    double x = cos(phi_val)*r_cone;
                    double y = sin(phi_val)*r_cone;
                    double _z = -r_cone*(_r_top-this->_r_cor)/_h_cor;

                    r_cone_arr.push_back(r_cone); da_arr.push_back(da);
                    x_arr.push_back(x); y_arr.push_back(y); z_arr.push_back(_z);

                    double ypos = -y;
                    double zpos = -z_val;
                    
                    ypos_arr.push_back(ypos); zpos_arr.push_back(zpos);

                    phi_val += dphi;
                }
                z_val += dz;
            }

            typename T<U>::iterator r_cone, da, x, y, z, ypos, zpos;

            for (omega = this->_omega_cor.begin(),
                    fcd = this->_frad_cortodisk.begin(),
                    fdc = this->_frad_disktocor.begin(), r = this->_r.begin();
                    omega != this->_omega_cor.end(),
                    fcd != this->_frad_cortodisk.end(),
                    fdc != this->_frad_disktocor.end(), r != this->_r.end();
                    omega++, fcd++, fdc++, r++) {
                double omegawt = 0;

                for (r_cone = r_cone_arr.begin(), da = da_arr.begin(),
                        x = x_arr.begin(), y = y_arr.begin(), z = z_arr.begin(),
                        ypos = ypos_arr.begin(), zpos = zpos_arr.begin();
                        r_cone != r_cone_arr.end(), da != da_arr.end(),
                        x != x_arr.end(), y != y_arr.end(), z != z_arr.end(),
                        ypos != ypos_arr.end(), zpos != zpos_arr.end();
                        r_cone++, da++, x++, y++, z++, ypos++, zpos++) {

                    double xpos = (*r)-(*x);
                    double rpos = pow(pow(xpos, 2) + pow(*ypos,2) +
                                        pow(*zpos, 2), .5);

                    double rcorpos = pow(pow(*r_cone, 2) + pow(*z, 2), .5);
                    double daomega = (*da)*((*x)*(xpos/rpos) +
                        (*y)*(*ypos/rpos) + (*z)*(*zpos/rpos))/
                        (pow(rpos, 2)*rcorpos);
                    double daomegawt = daomega*abs(*zpos/rpos);

                    if (daomega > 0 && *r > this->_r_cor) {
                        *omega += daomega;
                    }
                    if (daomegawt > 0 && *r > this->_r_cor) {
                        omegawt += daomegawt;
                    }
                }

                *fdc = omegawt/M_PI;
                *fcd = *fdc/(M_PI*(pow(_r_top, 2) + (this->_r_cor + _r_top) *
                    pow(pow((_r_top - this->_r_cor),2)+pow(_h_cor, 2), .5)));
            }
        }
};
