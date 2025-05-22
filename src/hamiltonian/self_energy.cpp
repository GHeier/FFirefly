#include "../config/load/cpp_config.hpp"
#include "../config/load/jl_interface.h"
#include "self_energy.hpp"

extern "C" void self_energy_wrapper() {
    call_self_energy();
}

void call_self_energy() {
    string folder = "hamiltonian/";
    string filename = "self_energy";
    string module = "Self_Energy";
    string function = "get_self_energy";
    call_julia_func(folder.c_str(), filename.c_str(), module.c_str(), function.c_str());
}

