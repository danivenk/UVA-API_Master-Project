
#include "../includes.h"

#ifndef geometry_H
#define geometry_H

template <class T, class U>
class Geometry {

    public:
        // constructor & destructor
        Geometry(T r, T r_area, U parms) : _r(r), _r_area(r_area),
                _parms(parms) {
            init(r, parms);
            setup_geometry();
            initialize_area_lists();
        };
        ~Geometry() {};
    
    private:
        // initialization
        void init(T r, U parms) {
            assert(tuple_size<U>{} >= 1);

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
        T _r;
        T _r_area;
        T _omega_cor, _frad_disktocor, _frad_cortodisk;
        T _omega_cor_area, _frad_disktocor_area, _frad_cortodisk_area;
        U _parms;
        
        void initialize_area_lists() {
            typename T::iterator omega, omega_area, fdc, fdc_area, fcd,
                fcd_area, area;

            for (omega = _omega_cor.begin(),
                    omega_area = _omega_cor_area.begin(),
                    fdc = _frad_disktocor.begin(),
                    fdc_area = _frad_disktocor_area.begin(),
                    fcd = _frad_cortodisk.begin(),
                    fcd_area = _frad_cortodisk_area.begin(),
                    area = _r_area.begin(); omega != _omega_cor.end(),
                    omega_area != _omega_cor_area.end(),
                    fdc != _frad_disktocor.end(),
                    fdc_area != _frad_disktocor_area.end(),
                    fcd != _frad_cortodisk.end(),
                    fcd_area != _frad_cortodisk_area.end(),
                    area != _r_area.end(); omega++, omega_area++, fdc++,
                    fdc_area++, fcd++, fcd_area++, area++) {
                *omega_area = (*omega)*(*area);
                *fdc_area = (*fdc)*(*area);
                *fcd_area = (*fcd)*(*area);
            }
        }

    public:
        // getters
        T get_rcor() { return get<0>(_parms); };
        T get_omega_cor() { return _omega_cor; };
        T get_frad_disktocor() { return _frad_disktocor; };
        T get_frad_cortodisk() { return _frad_cortodisk; };
        T get_omega_cor_area() { return _omega_cor_area; };
        T get_frad_disktocor_area() { return _frad_disktocor_area; };
        T get_frad_cortodisk_area() { return _frad_cortodisk_area; };
};

#endif