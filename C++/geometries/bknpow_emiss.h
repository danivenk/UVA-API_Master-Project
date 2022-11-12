#include "geometry.h"

template <template <typename...> class T, typename U>
class BKNpow_Emiss : public Geometry<T, U> {
    public:
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
        double total_frac = 0;
        U _r_bk1, _r_bk2, _em_ind1, _em_ind2, _em_ind3, _tot_ctd;

        void setup_geometry() {
            
            typename T<U>::iterator omega, fcd, fdc, r, area;

            for (omega = this->_omega_cor.begin(),
                    fcd = this->_frad_cortodisk.begin(),
                    fdc = this->_frad_disktocor.begin(), r = this->_r.begin(),
                    area = this->_r_area.begin();
                    omega != this->_omega_cor.end(),
                    fcd != this->_frad_cortodisk.end(),
                    fdc != this->_frad_disktocor.end(), r != this->_r.end(),
                    area != this->_r_area.end(); omega++, fcd++, fdc++, r++,
                    area++) {
                if (*r < _r_bk1 && *r > this->_r_cor) {
                    *fcd = pow((*r/_r_bk1), _em_ind1);
                } else if (*r >= _r_bk1 && *r < _r_bk2 && *r > this->_r_cor) {
                    *fcd = pow((*r/_r_bk1), _em_ind2);
                } else if (*r >= _r_bk2 && *r > this->_r_cor) {
                    *fcd = pow((_r_bk2/_r_bk1), _em_ind2) * 
                        pow((*r/_r_bk2), _em_ind3);
                }

                *fdc = *fcd*2*M_PI*pow(this->_r_cor, 2);

                total_frac += (*fcd)*(*area);
            }
        }

        void normalize_geometry() {
            typename T<U>::iterator fcd, fdc;

            for (fcd = this->_frad_cortodisk_area.begin(),
                    fdc = this->_frad_disktocor.begin();
                    fcd != this->_frad_cortodisk_area.end(),
                    fdc != this->_frad_disktocor.end();
                    fcd++, fdc++) {
                *fcd = (*fcd) * _tot_ctd / total_frac;
                *fdc = (*fdc) * _tot_ctd / total_frac;
            }
        }
};