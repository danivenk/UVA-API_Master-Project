// include the global includes
#include "../includes.h"

// make sure it is only included once
#ifndef geometry_H
#define geometry_H

/**
 * @brief Describes a geometry which can be used in the calc_illumination_frac
 *        function.
 *        It calculates the solid angle, corona to disk fraction & disk to
 *        corona fraction for each radial bin.
 * 
 * @tparam T is the container type
 * @tparam U is the value type
 */
template <template <typename...> class T, typename U>
class Geometry {

    public:
        /**
         * @brief Constructor for a Geometry class
         * 
         * @param r radial bins;
         * @param r_area radial area bins;
         * @param r_cor radius of the corona;
         */
        Geometry(T<U> r, T<U> r_area, U r_cor) : _r(r), _r_area(r_area),
                _r_cor(r_cor) {
            init(r);
            setup_geometry();
            initialize_area_lists();
        };
        ~Geometry() {};
    
    private:
        /** @brief Initialization of the class, setting up the arrays initially
         *         with only zeros.
         *  
         *  @param r radial bins;
         */
        void init(T<U> r) {

            // for each radial bin set the value of the arrays to 0
            for (int i = 0; i < r.size(); i++) {
                _omega_cor.push_back(0);
                _frad_cortodisk.push_back(0);
                _frad_disktocor.push_back(0);
                _omega_cor_area.push_back(0);
                _frad_cortodisk_area.push_back(0);
                _frad_disktocor_area.push_back(0);
            }
        };

        /**
         * @brief Set up the geometry
         */
        virtual void setup_geometry() {};

    protected:
        T<U> _r; // radial bins;
        T<U> _r_area; // radial area bins;
        // solid angle, disk to corona fraction, corona to disk fraction;
        T<U> _omega_cor, _frad_disktocor, _frad_cortodisk;
        T<U> _omega_cor_area, _frad_disktocor_area, _frad_cortodisk_area;
        U _r_cor; // corona radius
        
        /**
         * @brief Initialize the area lists
         */
        void initialize_area_lists() {

            // create the iterators needed
            typename T<U>::iterator omega, omega_area, fdc, fdc_area, fcd,
                fcd_area, area;

            // loop over arrays that need to be initialized 
            for (omega = _omega_cor.begin(),
                    omega_area = _omega_cor_area.begin(),
                    fdc = _frad_disktocor.begin(),
                    fdc_area = _frad_disktocor_area.begin(),
                    fcd = _frad_cortodisk.begin(),
                    fcd_area = _frad_cortodisk_area.begin(),
                    area = _r_area.begin(); omega != _omega_cor.end(),
                    omega_area != _omega_cor_area.end(); omega++, omega_area++,
                    fdc++, fdc_area++, fcd++, fcd_area++, area++) {
                // initialize values
                *omega_area = (*omega)*(*area);
                *fdc_area = (*fdc)*(*area);
                *fcd_area = (*fcd)*(*area);
            }
        }

    public:
        /**
         * @brief Get the radius of the corona
         * 
         * @return radius of the corona
         */
        U get_rcor() { return _r_cor; };

        /**
         * @brief Get the solid angle array
         * 
         * @return solid angle array
         */
        T<U> get_omega_cor() { return _omega_cor; };

        /**
         * @brief Get the disk to corona fraction array
         * 
         * @return disk to corona fraction array
         */
        T<U> get_frad_disktocor() { return _frad_disktocor; };

        /**
         * @brief Get the corona to disk fraction array
         * 
         * @return corona to disk fraction array
         */
        T<U> get_frad_cortodisk() { return _frad_cortodisk; };

        /**
         * @brief Get the solid angle times area array
         * 
         * @return solid angle times area array
         */
        T<U> get_omega_cor_area() { return _omega_cor_area; };

        /**
         * @brief Get the disk to corona fraction times area array
         * 
         * @return disk to corona fraction times area array
         */
        T<U> get_frad_disktocor_area() { return _frad_disktocor_area; };

        /**
         * @brief Get the corona to disk fraction times array
         * 
         * @return corona to disk fraction times array
         */
        T<U> get_frad_cortodisk_area() { return _frad_cortodisk_area; };
};

#endif