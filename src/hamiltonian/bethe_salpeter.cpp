#include "../config/load/cpp_config.hpp"
#include "../config/load/jl_interface.h"
#include "bethe_salpeter.hpp"

extern "C" void bethe_salpeter_wrapper() {
    call_BSE();
}

void call_BSE() {
    string folder = "hamiltonian/";
    string filename = "BSE";
    string module = "BSE";
    string function = "BSE_node";
    call_julia_func(folder.c_str(), filename.c_str(), module.c_str(), function.c_str());
}

