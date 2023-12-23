#include "functions.h"
#include <fstream>
#include <json.hpp>
#include <ctime>
#include <cstdlib>
#include <sstream>

using json = nlohmann::json;
// namespace py = pybind11;

template <class T>
ostream& operator<<(ostream& os, const list<T> list);
template <class T>
ostream& operator<<(ostream& os, const Matrix<list, T> m);
int myrandom (int i);
template <class T>
json to_json(complex<T> c);
template <class T>
list<json> to_json(list<complex<T>> c_list);
template <class T>
string to_string(complex<T> c);
template <class T>
list<string> to_string(list<complex<T>> c_list);

int main() {
    // time(0)
    srand (unsigned (0));
    
    list<double> array_1, array_2;
    double value_1 = .37;
    double value_2 = 6.7;

    for (int i = 1; i <= 10; i++) {
        array_1.push_back((double) i/10);
        array_2.push_front((double) i/10);
    }

    auto [dispfrac, seed_frac_flow, heat_frac_flow] =
        calc_dispfrac<double>(array_1, array_2, value_1, value_1, value_2,
                              value_1, value_1, value_1);

    return 0;
}

int myrandom (int i) { return std::rand()%i;}

template <class T>
json to_json(complex<T> c) {
    return {{"real", c.real()}, {"imag", c.imag()}};
}

template <class T>
list<json> to_json(list<complex<T>> c_list) {
    typename list<complex<T>>::iterator c;
    list<json> out;

    for (c = c_list.begin(); c != c_list.end(); c++) {
        out.push_back(to_json<T>(*c));
    }

    assert(c_list.size() == out.size());

    return out;
}

template <class T>
string to_string(complex<T> c) {
    stringstream ss;
    ss << c.real();
    if (c.imag() >= 0) {
        ss << "+";
    }
    ss << c.imag() << "j";
    return ss.str();
}

template <class T>
list<string> to_string(list<complex<T>> c_list) {
    typename list<complex<T>>::iterator c;
    list<string> out;
    
    for (c = c_list.begin(); c != c_list.end(); c++) {
        out.push_back(to_string<T>(*c));
    }
    
    assert(c_list.size() == out.size());
    
    return out;
}

template <class T>
ostream& operator<<(ostream& os, const list<T> list) {
    for (auto& i: list) {
        os << i << " ";
    }

    return os;
}

template <class T>
ostream& operator<<(ostream& os, const Matrix<list, T> m) {
    auto [x, y] = m.get_size();

    for (int i = 0; i < x; i++) {
        for (int j = 0; j < y; j++) {
            os << m.get_element(i, j) << " ";
        }
        os << endl;
    }

    return os;
}