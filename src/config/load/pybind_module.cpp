#include <pybind11/pybind11.h>
#include <pybind11/stl.h>  // For STL containers like std::string
#include "../../hamiltonian/potential.hpp"
#include "../../hamiltonian/band_structure.hpp"
#include "../../objects/vec.hpp"

namespace py = pybind11;

PYBIND11_MODULE(fcode, m) {
    // Expose the Vec class
    py::class_<Vec>(m, "Vec")
        .def(py::init<float, float, float>(), py::arg("x"), py::arg("y"), py::arg("z"))
        .def("norm", &Vec::norm, "Calculate the norm of the vector")
        .def("__getitem__", [](const Vec& v, size_t i) -> float {
            if (i == 0) return v.x;
            else if (i == 1) return v.y;
            else if (i == 2) return v.z;
            else throw py::index_error("Index out of range");
        })  // Allow indexing with []
        .def("__repr__", [](const Vec& v) {
            return "<Vec x=" + std::to_string(v.x) +
                   " y=" + std::to_string(v.y) +
                   " z=" + std::to_string(v.z) + ">";
        });

    // Expose the V function
    m.def("V", &V, "Calculate the interaction potential",
          py::arg("k1"), py::arg("k2"), py::arg("spin1"), py::arg("spin2"));

    // Expose the epsilon function
    m.def("epsilon", &epsilon, "Calculate the energy band",
          py::arg("n"), py::arg("k"));
}
