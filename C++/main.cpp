#include <iostream>
#include <list>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <tuple>
#include <cassert>
#include <iterator>
#include <typeinfo>

using namespace std;


template <class T>
tuple<T, int> find_nearest(list<T> array, T value);
template <class T>
ostream& operator<<(ostream& os, const list<T> list);

template <class T, class U>
class Geometry {

    public:
        // constructor & destructor
        Geometry(T r, list<U> parms) : _r(r), _parms(parms) {
            init(r, parms);
            setup_geometry();
        };
        ~Geometry() {};
    
    private:
        // initialization
        void init(T r, list<U> parms) {
            assert(parms.size() >= 1);

            for (int i = 0; i < r.size(); i++) {
                _omega_cor.push_back(0);
                _frad_cortodisk.push_back(0);
                _frad_disktocor.push_back(0);
            }
        };
        virtual void setup_geometry() {};

    protected:
        // parameter declaration
        T _r;
        T _omega_cor, _frad_disktocor, _frad_cortodisk;
        list<U> _parms;

    public:
        // getters
        T get_omega_cor() { return _omega_cor; };
        T get_frad_disktocor() { return _frad_disktocor; };
        T get_frad_cortodisk() { return _frad_cortodisk; };
};

template <class T, class U>
class BKNpow_Emiss : public Geometry<T, U> {
    public:
        BKNpow_Emiss(T r, list<U> parms) : Geometry<T, U>(r, parms) {};

    private:
        void setup_geometry() {
            for (int i = 0; i < _r.size(); i++) {
                _omega_cor.push_back(1);
                _frad_cortodisk.push_back(1);
                _frad_disktocor.push_back(1);
            }
        }
};

int main() {

    srand(std::time(nullptr));

    list<double> _list = {};

    for (int i = 0; i < 10; i++) {
        _list.push_back(rand());
    }

    Geometry<list<double>, double> a(_list, {1,2,3});

    cout << a.get_omega_cor() << endl;
    cout << a.get_frad_disktocor() << endl;
    cout << a.get_frad_cortodisk() << endl;
    
    BKNpow_Emiss<list<double>, double> b(_list, {1,2,3});

    cout << b.get_omega_cor() << endl;
    cout << b.get_frad_disktocor() << endl;
    cout << b.get_frad_cortodisk() << endl;

    // A a(2);

    // cout << "A" << endl;
    // cout << a.get_a() << endl;
    // cout << square(a) << endl;

    // B b(2);

    // cout << "B" << endl;
    // cout << b.get_a() << endl;
    // cout << square(b) << endl;

    return 0;
}

template <class T>
tuple<T, int> find_nearest(list<T> array, T value) {
    
    cout << "list: " << array << endl;
    cout << "value: " << value << endl;

    T smallest = numeric_limits<T>::max();
    T small_diff = numeric_limits<T>::max();

    int i = 0;

    for (T &item: array) {
        T diff = abs(item - value);

        cout << smallest << " - " << item << endl;
        cout << small_diff << " _ " << diff << endl;

        if (diff < small_diff) {
            smallest = item;
            small_diff = diff;
        }

        i++;
    }
    
    return make_tuple(smallest, i);
}

template <class T>
tuple<list<T>, list<T>, list<T>> calc_dispfrac(list<T> rad, list<T> rad_area,
        T rin, T rcor, T seedff_norm, T seedff_ind, T heatff_norm,
        T heatff_ind) {

    assert (rad.size() == rad_area.size());

    list<T> dispfrac (rad.size(), 0);
    list<T> seed_frac_flow = (rad.size(), 0);
    list<T> heat_frac_flow = (rad.size(), 0);

    typename list<T>::iterator disp, seed, heat, r, area;

    for (disp = dispfrac.begin(), seed = seed_frac_flow.begin(),
            heat = heat_frac_flow.begin(), r = rad.begin(),
            area = rad_area.begin(); disp != dispfrac.end(); disp++, seed++,
            heat++, r++, area++) {
        *disp = *area * pow(*r, -3.0) * (1 - sqrt(rin/(*r)));

        if (*r <= rcor) {
            *heat = pow(heatff_norm*(*r/rin), heatff_ind);

            if (seedff_norm >= 0) {
                *seed = pow(seedff_norm*(*r/rin), heatff_ind);
            } else {
                *seed = 1.0 - *heat;
            }
        }
    }

    return make_tuple(dispfrac, seed_frac_flow, heat_frac_flow);
}

template <class T>
tuple<list<T>, list<T>, list<T>> calc_illumination_fracs(list<T> rad,
        list<T> rad_area, Geometry<list<T>, T> geomod, list<T> parms) {

    assert(rad.size() == rad_area.size());

    list<T> omega_cor(rad.size(), 0), frad_disktocor(rad.size(), 0),
        frad_cortodisk(rad.size(), 0);



    return make_tuple(omega_cor, frad_disktocor, frad_cortodisk);
}

template <class T>
ostream& operator<<(ostream& os, const list<T> list) {
    for (auto& i: list) {
        os << i << " ";
    }

    return os;
}