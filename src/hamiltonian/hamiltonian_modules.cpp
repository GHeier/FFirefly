#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "potential.h"
#include "band_structure.h"

namespace py = pybind11;

using namespace std;

// Expose functions and types to Python
PYBIND11_MODULE(fcode, m) {
    py::class_<Vec>(m, "Vec")
        .def(py::init<float, float, float>())
        .def_readwrite("x", &Vec::x)
        .def_readwrite("y", &Vec::y)
        .def_readwrite("z", &Vec::z);

    m.def("epsilon", &epsilon, "epsilon(n,k)",
          py::arg("n"), py::arg("k"));

    m.def("V", &V, "V(k1, k2, spin1, spin2)",
          py::arg("k1"), py::arg("k2"), py::arg("spin1"), py::arg("spin2"));
}
