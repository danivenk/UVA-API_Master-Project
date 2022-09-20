#include "geometry.h"

template <class T, class U>
class Inv_Cone : public Geometry<T, U> {
    public:
        Inv_Cone(T r, T r_area, U parms) : Geometry<T, U>(r, r_area, parms) {
            setup_geometry();
            this->initialize_area_lists();
        };

    private:
        void setup_geometry() {
            
            typename T::iterator omega, fcd, fdc, r;

            assert(tuple_size<U>{} == 5);

            auto [r_cor, h_cor, r_top, nz, nphi] = this->_parms;

            double dz = h_cor/nz, dphi = 2*M_PI/nphi;

            T r_cone_arr, da_arr, x_arr, y_arr, z_arr,
                ypos_arr, zpos_arr;
                
            double cone_angle = atan((r_top-r_cor)/h_cor);

            for (double z = 0; z/dz < nz; z += dz) {
                for (double phi = 0; phi/dphi < nphi; phi += dphi) {
                    double r_cone = r_cor + z*(r_top-r_cor)/h_cor;
                    double da = dphi*r_cone*dz/cos(cone_angle);

                    double x = cos(phi)*r_cone;
                    double y = sin(phi)*r_cone;
                    double _z = -r_cone*(r_top-r_cor)/h_cor;

                    r_cone_arr.push_back(r_cone); da_arr.push_back(da);
                    x_arr.push_back(x); y_arr.push_back(y); z_arr.push_back(_z);

                    double ypos = -y;
                    double zpos = -z;
                    
                    ypos_arr.push_back(ypos); zpos_arr.push_back(zpos);
                }
            }

            typename T::iterator r_cone, da, x, y, z, ypos, zpos;

            for (omega = this->_omega_cor.begin(),
                    fcd = this->_frad_cortodisk.begin(),
                    fdc = this->_frad_disktocor.begin(), r = this->_r.begin();
                    omega != this->_omega_cor.end(),
                    fcd != this->_frad_cortodisk.end(),
                    fdc != this->_frad_disktocor.end(), r != this->_r.end();
                    omega++, fcd++, fdc++, r++) {
                double omegawt = 0;
                int i = 0;

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

                    if (daomega > 0 && *r > r_cor) {
                        *omega += daomega;
                    }
                    if (daomegawt > 0 && *r > r_cor) {
                        omegawt += daomegawt;
                    }
                }

                *fdc = omegawt/M_PI;
                *fcd = *fdc/(M_PI*(pow(r_top, 2) + (r_cor + r_top) *
                    pow(pow((r_top - r_cor),2)+pow(h_cor, 2), .5)));
            }
        }
};
