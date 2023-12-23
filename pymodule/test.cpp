#include <pybind11/embed.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>
#include <list>
#include "../C++/functions.h"
namespace py = pybind11;
using namespace std;

template<class T>
T myadd(T a, T b) {
    return a + b;
} 

PYBIND11_MODULE(example, m) {
    m.doc() = "TEST FUNCTION";
    m.def("myadd", &myadd<int>, "Add ints");
    m.def("myadd", &myadd<double>, "Add doubles");
    m.def("myadd", &myadd<std::string>, "Add strings");
    m.def("find_nearest", &find_nearest<double>);
}