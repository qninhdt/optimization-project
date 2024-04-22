// cppimport
#include <pybind11/pybind11.h>
#include "lmao.hpp"

namespace py = pybind11;

PYBIND11_MODULE(hello, m)
{
    m.def("square", &square);
}
/*

*/