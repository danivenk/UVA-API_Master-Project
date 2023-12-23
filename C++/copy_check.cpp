#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/embed.h>
#include <random>
#include <memory>
#include "includes.h"

namespace py = pybind11;
using namespace py::literals;

template <typename T, template <typename...> class U>
ostream& operator<<(ostream& os, const U<T> array) {
    for (auto& i: array) {
        os << i << " ";
    }

    return os;
}

template <typename T>
vector<T> vfa(py::array_t<T> array) {
    auto buf = array.request();
    T* ptr = static_cast<T*>(buf.ptr);
    return vector<T>(ptr, ptr + buf.size);
}

template <typename T>
py::array_t<T> afv(const std::vector<T>& vec) {
    py::array_t<T> arr(vec.size());
    T* ptr = arr.mutable_data();
    std::memcpy(ptr, vec.data(), vec.size() * sizeof(T));
    return arr;
}

template <typename T>
T myrandom (T i) { return std::rand()%i; }

int main() {
    int N = 100;
    double min(0), max(1);
    mt19937 rng(random_device{}());
    uniform_real_distribution<double> dist(min, max);
    vector<double> myVector;
    myVector.reserve(N);
    generate_n(back_inserter(myVector), N, [&]() { return dist(rng); });

    auto a = afv(myVector);

    return 0;
}