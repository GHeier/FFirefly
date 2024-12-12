#include "cpp_config.hpp"  
#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // For STL containers like std::string and arrays
#include "../../hamiltonian/potential.hpp"
#include "../../hamiltonian/band_structure.hpp"
#include "../../objects/vec.hpp"

namespace py = pybind11;

PYBIND11_MODULE(fcode, m) {

    // Bind the update function
    m.def("load_config", &load_cpp_config, 
          "Update the global variables");

    py::class_<Vec>(m, "Vec")
        .def(py::init<float, float, float>())
        .def_readwrite("x", &Vec::x)
        .def_readwrite("y", &Vec::y)
        .def_readwrite("z", &Vec::z);

    m.def("epsilon", &epsilon, "epsilon(n,k)",
          py::arg("n"), py::arg("k"));

    m.def("V", &V, "V(k1, k2, spin1, spin2)",
          py::arg("k1"), py::arg("k2"), py::arg("spin1"), py::arg("spin2"));


    // Begin PyBind Definitions

//[CONTROL]
    m.attr("category") = &category;
    m.attr("calculation") = &calculation;
    m.attr("outdir") = &outdir;
    m.attr("prefix") = &prefix;
    m.attr("verbosity") = &verbosity;
    m.attr("datfile_in") = &datfile_in;
    m.attr("datfile_out") = &datfile_out;

//[SYSTEM]
    m.attr("interaction") = &interaction;
    m.attr("dimension") = &dimension;
    m.attr("ibrav") = &ibrav;
    m.attr("nbnd") = &nbnd;
    m.attr("fermi_energy") = &fermi_energy;
    m.attr("Temperature") = &Temperature;
    m.attr("onsite_U") = &onsite_U;

//[MESH]
    m.attr("k_mesh") = &k_mesh;
    m.attr("q_mesh") = &q_mesh;
    m.attr("w_pts") = &w_pts;

//[CELL]
    m.attr("cell") = &cell;
    m.attr("brillouin_zone") = &brillouin_zone;

//[BANDS]
    m.attr("band") = &band;
    m.attr("eff_mass") = &eff_mass;
    m.attr("t0") = &t0;
    m.attr("t1") = &t1;
    m.attr("t2") = &t2;
    m.attr("t3") = &t3;
    m.attr("t4") = &t4;
    m.attr("t5") = &t5;
    m.attr("t6") = &t6;
    m.attr("t7") = &t7;
    m.attr("t8") = &t8;
    m.attr("t9") = &t9;
    m.attr("t10") = &t10;

//[SUPERCONDUCTOR]
    m.attr("FS_only") = &FS_only;
    m.attr("bcs_cutoff_frequency") = &bcs_cutoff_frequency;
    m.attr("num_eigenvalues_to_save") = &num_eigenvalues_to_save;
    m.attr("frequency_pts") = &frequency_pts;

//[RESPONSE]
    m.attr("dynamic") = &dynamic;
    // End PyBind Definitions

}
