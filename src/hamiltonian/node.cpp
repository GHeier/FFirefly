#include "node.hpp"
#include "../config/load/py_interface.h"
#include "../config/load/cpp_config.hpp"
#include "constant_energy_integration_DOS.hpp"
#include "fs.hpp"
#include <iostream>
#include <string>

using namespace std;

extern "C" void hamiltonian_wrapper() {
    printv("Running hamiltonian_wrapper\n");

    if (calculation == "generate") {
        generate_Hk();
    }
    else if (calculation == "DOS") {
        get_DOS();
    }
    else if (calculation == "FS") {
        save_FS();
    }
    else {
        printv("Hamiltonian category: calculation `%s` not recognized\n", calculation.c_str());
    }
}

void generate_Hk() {
    string folder = "hamiltonian/";
    string filename = "Hk";
    string function = "generate_hamiltonian";
    call_python_func(folder.c_str(), filename.c_str(), function.c_str());
}

void get_DOS() {
    if (method == "gaussian") {
        string folder = "hamiltonian/";
        string filename = "DOS";
        string function = "get_DOS";
        call_python_func(folder.c_str(), filename.c_str(), function.c_str());
    }
    else if (method == "tetrahedra") {
        tetrahedra_sum();
    }
    else {
        printv("DOS category: method `%s` not recognized\n", method.c_str());
    }
}
