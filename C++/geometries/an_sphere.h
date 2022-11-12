#include "geometry.h"

template <template <typename...> class T, typename U>
class AN_Sphere : public Geometry<T, U> {
    public:
        AN_Sphere(T<U> r, T<U> r_area, U r_cor) :
                Geometry<T, U>(r, r_area, r_cor) {
            setup_geometry();
            this->initialize_area_lists();
        };

    private:
        void setup_geometry() {
            
            typename T<U>::iterator fcd, fdc, r;

            for (fcd = this->_frad_cortodisk.begin(),
                    fdc = this->_frad_disktocor.begin(), r = this->_r.begin();
                    fcd != this->_frad_cortodisk.end(),
                    fdc != this->_frad_disktocor.end(), r != this->_r.end();
                    fcd++, fdc++, r++) {
                if (*r > this->_r_cor) {
                    *fdc = (asin(this->_r_cor/(*r)) - 
                        sqrt(pow((this->_r_cor/(*r)), 2) -
                        pow((this->_r_cor/(*r)), 4)))/M_PI;
                    *fcd = *fdc/(2*M_PI*pow(this->_r_cor, 2));
                }
            }
        }
};