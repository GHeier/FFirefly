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

float add(float i, float j) {
    return i + j;
}

float e_test(Vec v) {
    //return epsilon(0, v);
    return add(v.x, v.y);
}

extern "C" {
    void load_c_config();
}

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
        .def("norm", &Vec::norm);

    m.def("add", &add, "A function which adds two numbers");

    m.def("e_test", &e_test, "A function which returns the x component of a Vec");

    // Expose the V function
    m.def("V", &V, "Calculate the interaction potential",
          py::arg("k1"), py::arg("k2"), py::arg("spin1") = "up", py::arg("spin2") = "up");


    // Expose the epsilon function
    m.def("epsilon", &epsilon, "Calculate the energy band",
            py::arg("n") = 0, py::arg("k"));
    py::class_<Config>(m, "Config")
        .def(py::init<>())
    // Begin the Config class

//[CONTROL]
    .def_readwrite("category", &Config::category)
    .def_readwrite("calculation", &Config::calculation)
    .def_readwrite("outdir", &Config::outdir)
    .def_readwrite("prefix", &Config::prefix)
    .def_readwrite("verbosity", &Config::verbosity)
    .def_readwrite("datfile_in", &Config::datfile_in)
    .def_readwrite("datfile_out", &Config::datfile_out)

//[SYSTEM]
    .def_readwrite("interaction", &Config::interaction)
    .def_readwrite("dimension", &Config::dimension)
    .def_readwrite("ibrav", &Config::ibrav)
    .def_readwrite("nbnd", &Config::nbnd)
    .def_readwrite("fermi_energy", &Config::fermi_energy)
    .def_readwrite("Temperature", &Config::Temperature)
    .def_readwrite("onsite_U", &Config::onsite_U)

//[MESH]
    .def_readwrite("k_mesh", &Config::k_mesh)
    .def_readwrite("q_mesh", &Config::q_mesh)
    .def_readwrite("w_pts", &Config::w_pts)

//[CELL]
    .def_readwrite("cell", &Config::cell)

//[BRILLOUIN_ZONE]
    .def_readwrite("brillouin_zone", &Config::brillouin_zone)

//[BANDS]
    .def_readwrite("band", &Config::band)
    .def_readwrite("eff_mass", &Config::eff_mass)
    .def_readwrite("t0", &Config::t0)
    .def_readwrite("t1", &Config::t1)
    .def_readwrite("t2", &Config::t2)
    .def_readwrite("t3", &Config::t3)
    .def_readwrite("t4", &Config::t4)
    .def_readwrite("t5", &Config::t5)
    .def_readwrite("t6", &Config::t6)
    .def_readwrite("t7", &Config::t7)
    .def_readwrite("t8", &Config::t8)
    .def_readwrite("t9", &Config::t9)
    .def_readwrite("t10", &Config::t10)

//[SUPERCONDUCTOR]
    .def_readwrite("method", &Config::method)
    .def_readwrite("FS_only", &Config::FS_only)
    .def_readwrite("bcs_cutoff_frequency", &Config::bcs_cutoff_frequency)
    .def_readwrite("num_eigenvalues_to_save", &Config::num_eigenvalues_to_save)
    .def_readwrite("frequency_pts", &Config::frequency_pts)

//[RESPONSE]
    .def_readwrite("dynamic", &Config::dynamic)
    // End the Config class
    ;
    m.def("py_to_cpp_config", [](Config &c) {
        pc = &c;
    });
    m.def("load_c_config", &load_c_config);
    m.def("load_cpp_config", &load_cpp_config);
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
