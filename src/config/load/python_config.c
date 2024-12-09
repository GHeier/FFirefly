#include "c_config.h"  
#include <pybind11/pybind11.h>

PYBIND11_MODULE(fcode, m) {
    // Begin PyBind Definitions

//[CONTROL]
    m.attr("category") = &c_category;
    m.attr("calculation") = &c_calculation;
    m.attr("outdir") = &c_outdir;
    m.attr("prefix") = &c_prefix;
    m.attr("verbosity") = &c_verbosity;
    m.attr("datfile_in") = &c_datfile_in;
    m.attr("datfile_out") = &c_datfile_out;

//[SYSTEM]
    m.attr("interaction") = &c_interaction;
    m.attr("dimension") = &c_dimension;
    m.attr("ibrav") = &c_ibrav;
    m.attr("nbnd") = &c_nbnd;
    m.attr("fermi_energy") = &c_fermi_energy;
    m.attr("Temperature") = &c_Temperature;
    m.attr("onsite_U") = &c_onsite_U;

//[MESH]
    m.attr("k_mesh") = &c_k_mesh;
    m.attr("q_mesh") = &c_q_mesh;
    m.attr("w_pts") = &c_w_pts;

//[CELL]
    m.attr("cell") = &c_cell;
    m.attr("brillouin_zone") = &c_brillouin_zone;

//[BANDS]
    m.attr("band") = &c_band;
    m.attr("eff_mass") = &c_eff_mass;

//[SUPERCONDUCTOR]
    m.attr("FS_only") = &c_FS_only;
    m.attr("bcs_cutoff_frequency") = &c_bcs_cutoff_frequency;
    m.attr("num_eigenvalues_to_save") = &c_num_eigenvalues_to_save;
    m.attr("frequency_pts") = &c_frequency_pts;

//[RESPONSE]
    m.attr("dynamic") = &c_dynamic;
    // End PyBind Definitions

    // Bind the update function
    m.def("load_config", &load_c_config, 
          "Update the global variables");
}
