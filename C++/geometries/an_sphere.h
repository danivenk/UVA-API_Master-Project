#include "geometry.h"

template <class T, class U>
class AN_Sphere : public Geometry<T, U> {
    public:
        AN_Sphere(T r, T r_area, U parms) : Geometry<T, U>(r, r_area, parms) {
            setup_geometry();
            this->initialize_area_lists();
        };

    private:
        void setup_geometry() {
            
            typename T::iterator fcd, fdc, r;

            assert(tuple_size<U>{} == 1);

            auto [r_cor] = this->_parms;

            for (fcd = this->_frad_cortodisk.begin(),
                    fdc = this->_frad_disktocor.begin(), r = this->_r.begin();
                    fcd != this->_frad_cortodisk.end(),
                    fdc != this->_frad_disktocor.end(), r != this->_r.end();
                    fcd++, fdc++, r++) {
                if (*r > r_cor) {
                    *fdc = (asin(r_cor/(*r))-sqrt(pow((r_cor/(*r)), 2) -
                        pow((r_cor/(*r)), 4)))/M_PI;
                    *fcd = *fdc/(2*M_PI*pow(r_cor, 2));
                }
            }
        }
};