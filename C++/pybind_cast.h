#include <pybind11/pybind11.h>
// #include <pybind11/numpy.h>
// #include "functions.h"

namespace py = pybind11;

namespace pybind11 { namespace detail {

    template <>
    struct type_caster<int> {
    PYBIND11_TYPE_CASTER(int, _("int"));

    bool load(handle src, bool) {
        if (PyNumber_Check(src.ptr())) {
            auto value = PyLong_AsLong(src.ptr());
            return true;
        }
        return false;
    }

    static handle cast(int src, return_value_policy, handle) {
        return PyLong_FromLong(src);
    }
    };

}}