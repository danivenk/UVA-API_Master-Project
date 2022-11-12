
#include "../includes.h"

#ifndef geometry_H
#define geometry_H

template <template <typename...> class T, typename U>
class Geometry {

    public:
        // constructor & destructor
        Geometry(T<U> r, T<U> r_area, U r_cor) : _r(r), _r_area(r_area),
                _r_cor(r_cor) {
            init(r);
            setup_geometry();
            initialize_area_lists();
        };
        ~Geometry() {};
    
    private:
        // initialization
        void init(T<U> r) {

            for (int i = 0; i < r.size(); i++) {
                _omega_cor.push_back(0);
                _frad_cortodisk.push_back(0);
                _frad_disktocor.push_back(0);
                _omega_cor_area.push_back(0);
                _frad_cortodisk_area.push_back(0);
                _frad_disktocor_area.push_back(0);
            }
        };
        virtual void setup_geometry() {};

    protected:
        // parameter declaration
        T<U> _r;
        T<U> _r_area;
        T<U> _omega_cor, _frad_disktocor, _frad_cortodisk;
        T<U> _omega_cor_area, _frad_disktocor_area, _frad_cortodisk_area;
        U _r_cor;
        
        void initialize_area_lists() {
            typename T<U>::iterator omega, omega_area, fdc, fdc_area, fcd,
                fcd_area, area;

            for (omega = _omega_cor.begin(),
                    omega_area = _omega_cor_area.begin(),
                    fdc = _frad_disktocor.begin(),
                    fdc_area = _frad_disktocor_area.begin(),
                    fcd = _frad_cortodisk.begin(),
                    fcd_area = _frad_cortodisk_area.begin(),
                    area = _r_area.begin(); omega != _omega_cor.end(),
                    omega_area != _omega_cor_area.end(); omega++, omega_area++,
                    fdc++, fdc_area++, fcd++, fcd_area++, area++) {
                *omega_area = (*omega)*(*area);
                *fdc_area = (*fdc)*(*area);
                *fcd_area = (*fcd)*(*area);
            }
        }

    public:
        // getters
        U get_rcor() { return _r_cor; };
        T<U> get_omega_cor() { return _omega_cor; };
        T<U> get_frad_disktocor() { return _frad_disktocor; };
        T<U> get_frad_cortodisk() { return _frad_cortodisk; };
        T<U> get_omega_cor_area() { return _omega_cor_area; };
        T<U> get_frad_disktocor_area() { return _frad_disktocor_area; };
        T<U> get_frad_cortodisk_area() { return _frad_cortodisk_area; };
};

#endif