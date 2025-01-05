#include <pybind11/pybind11.h>
#include <pybind11/stl.h>  // For STL containers like std::string
                           //
#include <vector>
#include <string>


#include "config/load/cpp_config.hpp"
#include "config/load/c_config.h"
#include "hamiltonian/interaction.hpp"
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
        .def(py::init<vector<Vec>, vector<float>, int>())
        .def(py::init<vector<vector<double>>, vector<float>, int>())
        .def(py::init<string, int>())
        .def("__call__", py::overload_cast<vector<double>>(&ScalarField::operator()))
        .def("__call__", py::overload_cast<Vec>(&ScalarField::operator()));

    py::class_<ComplexField>(m, "ComplexField")
        .def(py::init<vector<Vec>, vector<complex<float>>, int>())
        .def(py::init<vector<vector<double>>, vector<complex<float>>, int>())
        .def(py::init<string, int>())
        .def("__call__", py::overload_cast<vector<double>>(&ComplexField::operator()))
        .def("__call__", py::overload_cast<Vec>(&ComplexField::operator()));

    py::class_<VectorField>(m, "VectorField")
        .def(py::init<vector<Vec>, vector<float>, int>())
        .def(py::init<vector<vector<double>>, vector<float>, int>())
        .def(py::init<string, int>())
        .def("__call__", py::overload_cast<vector<double>>(&VectorField::operator()))
        .def("__call__", py::overload_cast<Vec>(&VectorField::operator()));

    py::class_<ComplexVectorField>(m, "ComplexVectorField")
        .def(py::init<vector<Vec>, vector<complex<float>>, int>())
        .def(py::init<vector<vector<double>>, vector<complex<float>>, int>())
        .def(py::init<string, int>())
        .def("__call__", py::overload_cast<vector<double>>(&ComplexVectorField::operator()))
        .def("__call__", py::overload_cast<Vec>(&ComplexVectorField::operator()));

    // Expose the V function
    m.def("V", &V, "Calculate the interaction potential",
          py::arg("k1"), py::arg("k2"), py::arg("spin1") = "up", py::arg("spin2") = "up");


    // Expose the epsilon function
    m.def("epsilon", &epsilon, "Calculate the energy band",
            py::arg("n") = 1, py::arg("k"));

    m.def("epsilon", [](int n, py::array_t<double> k) -> py::object {
        if (k.ndim() == 1 && k.size() == 3) {
            // Single k-point case
            auto k_unchecked = k.unchecked<1>();
            Vec k_vec = Vec(k_unchecked(0), k_unchecked(1), k_unchecked(2));
            return py::float_(epsilon(n, k_vec));
        } 
        else if (k.ndim() == 2 && k.shape(1) == 3) {
            // Multiple k-points case
            auto buf = k.request();
            py::array_t<double> result(buf.shape[0]);
            auto r = result.mutable_unchecked<1>();
            auto k_unchecked = k.unchecked<2>();

            for (ssize_t i = 0; i < buf.shape[0]; i++) {
                Vec k_vec = Vec(k_unchecked(i, 0), k_unchecked(i, 1), k_unchecked(i, 2));
                r(i) = epsilon(n, k_vec);
            }
            return result;
        } 
        else {
            throw std::runtime_error("Input must have shape (3,) or (N, 3).");
        }
    });


    m.def("load_c_config", [](string path) {
        read_c_config_wrapper(path);
        load_cpp_config();
    });

}
