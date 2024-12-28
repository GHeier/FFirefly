#include <pybind11/pybind11.h>
#include <pybind11/stl.h>  // For STL containers like std::string
                           //
#include <vector>
#include <string>


#include "config/load/cpp_config.hpp"
#include "hamiltonian/potential.hpp"
#include "hamiltonian/band_structure.hpp"
#include "objects/vec.hpp"


namespace py = pybind11;

PYBIND11_MODULE(fmodule, m) {
    // Expose the Vec class
    py::class_<Vec>(m, "Vec")
        // Expose the default constructor
        .def(py::init<>())
        // Expose the parameterized constructor
        .def(py::init<float, float, float, float, float, int, int>(),
             py::arg("x"), py::arg("y") = 0.0f, py::arg("z") = 0.0f,
             py::arg("w") = 0.0f, py::arg("area") = 0.0f,
             py::arg("dimension") = 3, py::arg("n") = 1)
        // Optionally expose the class members as read/write attributes
        .def_readwrite("x", &Vec::x)
        .def_readwrite("y", &Vec::y)
        .def_readwrite("z", &Vec::z)
        .def_readwrite("w", &Vec::w)
        .def_readwrite("area", &Vec::area)
        .def_readwrite("dimension", &Vec::dimension)
        .def_readwrite("n", &Vec::n)
        // Norm function
        .def("norm", &Vec::norm)

        .def("__mul__", [](const Vec &v, int m) { return v * m; }, py::arg("multiple"))
        .def("__mul__", [](const Vec &v, float m) { return v * m; }, py::arg("multiple"))
        .def("__rmul__", [](int m, const Vec &v) { return m * v; }, py::arg("multiple"))
        .def("__rmul__", [](float m, const Vec &v) { return m * v; }, py::arg("multiple"))
        // Binding __truediv__ (division) for both int and float
        .def("__truediv__", [](const Vec &v, int d) { return v / d; }, py::arg("multiple"))
        .def("__truediv__", [](const Vec &v, float d) { return v / d; }, py::arg("multiple"))
        .def("__rtruediv__", [](int d, const Vec &v) { return d / v; }, py::arg("multiple"))
        .def("__rtruediv__", [](float d, const Vec &v) { return d / v; }, py::arg("multiple"))

        .def("__call__", [](Vec &self, int i) -> float& {
                return self(i);
            }, py::return_value_policy::reference, py::arg("i"));

    py::class_<ScalarField>(m, "ScalarField")
        .def(py::init<>())
        .def(py::init<std::string, int, bool>(),
             py::arg("filename"), py::arg("dimension") = 3, py::arg("is_complex") = false)
        .def(py::init<std::vector<Vec>, std::vector<float>, int, bool>(),
             py::arg("points"), py::arg("values"), py::arg("dimension") = 3, py::arg("is_complex") = false)
        .def_readwrite("points", &ScalarField::points)
        .def_readwrite("values", &ScalarField::values)
        .def("__call__", 
             static_cast<float (ScalarField::*)(Vec)>(&ScalarField::operator()), 
             "Evaluate the scalar field at a Vec point")
        .def("__call__", 
             static_cast<float (ScalarField::*)(py::array_t<double>)>(&ScalarField::operator()), 
             "Evaluate the scalar field at a NumPy array");        

    // Expose the V function
    m.def("V", &V, "Calculate the interaction potential",
          py::arg("k1"), py::arg("k2"), py::arg("spin1") = "up", py::arg("spin2") = "up");


    // Expose the epsilon function
    m.def("epsilon", &epsilon, "Calculate the energy band",
            py::arg("n") = 0, py::arg("k"));
    ;

    m.def("set_nbnd", [](int nbnd_) {
        nbnd = nbnd_;
    });
    m.def("get_nbnd", []() {
        return nbnd;
    });
    m.def("set_onsite_U", [](float onsite_U_) {
        onsite_U = onsite_U_;
    });
    m.def("get_onsite_U", []() {
        return onsite_U;
    });
    m.def("set_category", [](std::string category_) {
        category = category_;
    });
    m.def("set_FS_only", [](bool FS_only_) {
        FS_only = FS_only_;
    });
    m.def("set_k_mesh", [](std::vector<int> &k_mesh_) {
        k_mesh = k_mesh_;
    });
    m.def("get_k_mesh", []() {
        return k_mesh;
    });
}
