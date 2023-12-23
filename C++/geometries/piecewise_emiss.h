#include "geometry.h"

template <template <typename...> class T, typename U>
class Piecewise_Emiss : public Geometry<T, U> {
    public:
        template<typename... Args>
        Piecewise_Emiss(T<U> r, T<U> r_area, U r_cor, U r_outer, Args&&... args) :
                Geometry<T, U>(r, r_area, r_cor) {
            _r_outer = r_outer;
            assert(sizeof...(args) % 3 == 0); // && sizeof...(args) > 5);
            init_args(forward<Args>(args)...);
            setup_geometry();
            this->initialize_area_lists();
            normalize_geometry();
        };

    private:
        U _r_outer;
        T<U> _r_inners, _em_inds, _norms;
        double total_frac = 0;

        void init_args(){};
        template<typename T_, typename ...Args_>
        void init_args(T_&& arg, Args_&&... args){
            cout << arg << endl;
            switch ((_r_inners.size() + _em_inds.size() + _norms.size()) % 3) {
                case 0:
                    _r_inners.push_back(forward<T_>(arg)); break;
                case 1:
                    _em_inds.push_back(forward<T_>(arg)); break;
                case 2:
                    _norms.push_back(forward<T_>(arg)); break;
            }
            init_args(forward<Args_>(args)...);
        };

        void setup_geometry() {
            
            typename T<U>::iterator omega, fcd, fdc, r, area, r_inner, em_ind,
                norm;

            cout << "list " << _em_inds << endl;
            cout << "begin " << *_em_inds.begin() << endl;

            for (omega = this->_omega_cor.begin(),
                    fcd = this->_frad_cortodisk.begin(),
                    fdc = this->_frad_disktocor.begin(), r = this->_r.begin(),
                    area = this->_r_area.begin();
                    omega != this->_omega_cor.end(),
                    fcd != this->_frad_cortodisk.end(),
                    fdc != this->_frad_disktocor.end(), r != this->_r.end(),
                    area != this->_r_area.end(); omega++, fcd++, fdc++, r++,
                    area++) {
                for (r_inner = _r_inners.begin(), em_ind = _em_inds.begin(),
                        norm = _norms.begin(); r_inner != _r_inners.end(),
                        em_ind != _em_inds.end(), norm != _norms.end();
                        r_inner++, em_ind++, norm++) {
                    if ((*r < _r_outer) && (*r >= *r_inner) && *r > this->_r_cor) {
                        cout << "ih" << endl;
                        *fdc = *norm * pow(*r / *r_inner, *em_ind);
                        cout << *fdc << " " << *norm << " " << *r << " " << *r_inner << " " << *em_ind << endl;
                        cout << this->_frad_disktocor << endl;
                    }
                    _r_outer = *r_inner;
                }

                *fcd = *fdc/(2*M_PI*pow(this->_r_cor, 2));

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
                *fcd = (*fcd) * (*_r_inners.begin()) / total_frac;
            }
        }
};