#include "geometry.h"

template <class T, class U>
class BKNpow_Emiss : public Geometry<T, U> {
    public:
        BKNpow_Emiss(T r, T r_area, U parms) : Geometry<T, U>(r, r_area, parms) {
            setup_geometry();
            this->initialize_area_lists();
            normalize_geometry();
        };

    private:
        double total_frac = 0;

        void setup_geometry() {
            
            typename T::iterator omega, fcd, fdc, r, area;

            assert(tuple_size<U>{} == 7);

            auto [r_cor, r_bk1, r_bk2, em_ind1, em_ind2, em_ind3,
                  total_cortodisk] = this->_parms;

            for (omega = this->_omega_cor.begin(),
                    fcd = this->_frad_cortodisk.begin(),
                    fdc = this->_frad_disktocor.begin(), r = this->_r.begin(),
                    area = this->_r_area.begin();
                    omega != this->_omega_cor.end(),
                    fcd != this->_frad_cortodisk.end(),
                    fdc != this->_frad_disktocor.end(), r != this->_r.end(),
                    area != this->_r_area.end(); omega++, fcd++, fdc++, r++,
                    area++) {
                if (*r < r_bk1 && *r > r_cor) {
                    *fcd = pow((*r/r_bk1), em_ind1);
                } else if (*r >= r_bk1 && *r < r_bk2 && *r > r_cor) {
                    *fcd = pow((*r/r_bk1), em_ind2);
                } else if (*r >= r_bk2 && *r > r_cor) {
                    *fcd = pow((*r/r_bk2), em_ind2)*pow((*r/r_bk2), em_ind3);
                }

                *fdc = *fcd*2*M_PI*pow(r_cor, 2);

                total_frac += (*fcd)*(*area);
            }
        }

        void normalize_geometry() {
            typename T::iterator fcd, fdc;

            for (fcd = this->_frad_cortodisk_area.begin(),
                    fdc = this->_frad_disktocor_area.begin();
                    fcd != this->_frad_cortodisk_area.end(),
                    fdc != this->_frad_disktocor_area.end();
                    fcd++, fdc++) {
                *fcd = (*fcd) * get<6>(this->_parms) / total_frac;
                *fdc = (*fdc) * get<6>(this->_parms) / total_frac;
            }
        }
};