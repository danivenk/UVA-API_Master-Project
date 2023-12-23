#include <chrono>
#include <nlohmann/json.hpp>
#include <pybind11/pybind11.h>
#include <pybind11_json/pybind11_json.hpp>

namespace py = pybind11;
namespace nl = nlohmann;

using namespace pybind11::literals;

nl::json global_json;  // Global JSON object

void a() {
    auto start = std::chrono::high_resolution_clock::now();
    // Code for function a
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    global_json["a_time"] = duration.count();
}

void b() {
    auto start = std::chrono::high_resolution_clock::now();
    // Code for function b
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    global_json["b_time"] = duration.count();
}

PYBIND11_MODULE(your_module, m) {
    m.def("a", &a);
    m.def("b", &b);
    m.attr("global_json") = &global_json;
}